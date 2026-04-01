#pragma once
#include <set>
#include "grid.hpp"
#include "solver_result.hpp"

// Determine which catalyst cells are "relevant" to the catalytic interaction.
//
// Algorithm: forward-simulate the evolution treating all catalyst cells as
// "unknown but stable." Use mask-based constraint propagation to resolve
// as many cells as possible. When the CGOL outcome is still ambiguous,
// peek at the solution and mark those catalyst cells as relevant.
//
// SM namespace and StableMaskGrid are defined in grid.hpp (shared with
// the pre-processing forward simulation).

// ── Forward simulation with peeking ──

inline std::set<std::pair<int,int>> find_relevant_cells(
    const Grid& grid, const SolverResult& result)
{
    std::set<std::pair<int,int>> relevant;
    int T = grid.total_gens;
    int H = grid.height, W = grid.width;

    // Phase 1: mask-based stability propagation
    StableMaskGrid smg;
    smg.init(grid.catalyst_positions, grid.perturbation_region);
    smg.propagate();
    auto stable = smg.stable_map();

    // Peek function
    auto peek = [&](int wx, int wy) {
        auto it = stable.find({wx, wy});
        if (it == stable.end() || it->second != -1) return;
        relevant.insert({wx, wy});
        int var = grid.catalyst_var_at(wx, wy);
        bool alive;
        if (var >= 2) alive = result.solution.count(var - 1) > 0;
        else alive = (var == 1);
        it->second = alive ? 1 : 0;
    };

    // Phase 2: forward simulation
    // states[t][gy][gx]: 0=dead, 1=alive, -1=unknown
    bool restart;
    do {
        restart = false;

        std::vector<std::vector<std::vector<int>>> states(T,
            std::vector<std::vector<int>>(H, std::vector<int>(W, 0)));

        // Initialize t=0
        for (int gy = 0; gy < H; gy++) {
            for (int gx = 0; gx < W; gx++) {
                int val = grid.cells[0][gy][gx];
                if (val == 0) {
                    states[0][gy][gx] = 0;
                } else if (val == 1) {
                    states[0][gy][gx] = 1;
                } else {
                    int wx = gx + grid.ox, wy = gy + grid.oy;
                    auto it = stable.find({wx, wy});
                    states[0][gy][gx] = (it != stable.end()) ? it->second : 0;
                }
            }
        }

        for (int t = 0; t < T - 1 && !restart; t++) {
            for (int gy = 0; gy < H && !restart; gy++) {
                for (int gx = 0; gx < W && !restart; gx++) {
                    // Skip cells outside light cone at t+1
                    if (grid.cells[t+1][gy][gx] == 0) continue;

                    // Stability optimization: if this cell AND all 8 neighbors
                    // match the stable catalyst state (known or still unknown),
                    // the configuration is still the stable still-life and stays
                    // the same. Works even for partially-peeked neighborhoods:
                    // a peeked cell with value matching stable[] counts as stable,
                    // and an unknown cell (-1) is also stable (hasn't been perturbed).
                    {
                        int wx = gx + grid.ox, wy = gy + grid.oy;
                        bool in_pert = grid.perturbation_region.count({wx, wy}) > 0
                                    || grid.catalyst_positions.count({wx, wy}) > 0;
                        if (in_pert) {
                            auto matches_stable = [&](int ngx, int ngy, int val) -> bool {
                                if (val == -1) return true; // still unknown = unperturbed
                                int wgx = grid.wrap_x(ngx);
                                int wgy = grid.wrap_y(ngy);
                                int nwx = wgx + grid.ox, nwy = wgy + grid.oy;
                                auto it = stable.find({nwx, nwy});
                                if (it != stable.end()) {
                                    if (it->second == -1) return true; // stable unknown
                                    return val == it->second;
                                }
                                return val == 0;
                            };
                            bool all_stable = matches_stable(gx, gy, states[t][gy][gx]);
                            if (all_stable) {
                                for (int dy = -1; dy <= 1 && all_stable; dy++)
                                    for (int dx = -1; dx <= 1 && all_stable; dx++) {
                                        if (dx == 0 && dy == 0) continue;
                                        int ny = grid.wrap_y(gy+dy);
                                        int nx = grid.wrap_x(gx+dx);
                                        if (!matches_stable(nx, ny, states[t][ny][nx]))
                                            all_stable = false;
                                    }
                            }
                            if (all_stable) {
                                states[t+1][gy][gx] = states[t][gy][gx];
                                continue;
                            }
                        }
                    }

                    int center = states[t][gy][gx];
                    int alive_n = 0, unknown_n = 0;

                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            if (dx == 0 && dy == 0) continue;
                            int ny = grid.wrap_y(gy + dy);
                            int nx = grid.wrap_x(gx + dx);
                            int s = states[t][ny][nx];
                            if (s == 1) alive_n++;
                            else if (s == -1) unknown_n++;
                        }
                    }

                    // Fully determined
                    if (center != -1 && unknown_n == 0) {
                        if (center == 1)
                            states[t+1][gy][gx] = (alive_n == 2 || alive_n == 3) ? 1 : 0;
                        else
                            states[t+1][gy][gx] = (alive_n == 3) ? 1 : 0;
                        continue;
                    }

                    int lo = alive_n, hi = alive_n + unknown_n;
                    bool determined = false;
                    int next = 0;

                    if (center == 1) {
                        if (lo >= 2 && hi <= 3)    { next = 1; determined = true; }
                        else if (hi < 2 || lo > 3) { next = 0; determined = true; }
                    } else if (center == 0) {
                        if (lo == 3 && hi == 3)    { next = 1; determined = true; }
                        else if (lo > 3 || hi < 3) { next = 0; determined = true; }
                    } else {
                        // Unknown center
                        bool alive_all_survive = (lo >= 2 && hi <= 3);
                        bool alive_all_die     = (hi < 2 || lo > 3);
                        bool dead_all_born     = (lo == 3 && hi == 3);
                        bool dead_all_stay     = (lo > 3 || hi < 3);

                        if (alive_all_survive && dead_all_born) {
                            next = 1; determined = true;
                        } else if (alive_all_die && dead_all_stay) {
                            next = 0; determined = true;
                        } else if (alive_all_survive && dead_all_stay) {
                            // Result = center state. Still unknown.
                            states[t+1][gy][gx] = -1;
                            continue;
                        } else if (alive_all_die && dead_all_born) {
                            // Result = !center. Still unknown.
                            states[t+1][gy][gx] = -1;
                            continue;
                        }
                    }

                    if (determined) {
                        states[t+1][gy][gx] = next;
                        continue;
                    }

                    // Ambiguous — peek unknown cells
                    if (center == -1) {
                        int wx = gx + grid.ox, wy = gy + grid.oy;
                        if (grid.catalyst_positions.count({wx, wy})) {
                            peek(wx, wy);
                            restart = true;
                            continue;
                        }
                    }

                    for (int dy = -1; dy <= 1 && !restart; dy++) {
                        for (int dx = -1; dx <= 1 && !restart; dx++) {
                            if (dx == 0 && dy == 0) continue;
                            int ny = grid.wrap_y(gy + dy);
                            int nx = grid.wrap_x(gx + dx);
                            if (states[t][ny][nx] == -1) {
                                int wx = nx + grid.ox, wy = ny + grid.oy;
                                if (grid.catalyst_positions.count({wx, wy})) {
                                    peek(wx, wy);
                                    restart = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    } while (restart);

    return relevant;
}

// Generate a smart blocking clause: negate only relevant catalyst cells' values.
inline std::vector<int> smart_blocking_clause(
    const Grid& grid, const SolverResult& result)
{
    auto relevant = find_relevant_cells(grid, result);

    std::vector<int> clause;
    for (auto [wx, wy] : relevant) {
        int var = grid.catalyst_var_at(wx, wy);
        if (var < 2) continue;
        int sat_var = var - 1;
        bool alive = result.solution.count(sat_var) > 0;
        clause.push_back(alive ? -sat_var : sat_var);
    }

    if (clause.empty()) {
        for (auto [wx, wy] : grid.catalyst_positions) {
            int var = grid.catalyst_var_at(wx, wy);
            if (var < 2) continue;
            int sat_var = var - 1;
            bool alive = result.solution.count(sat_var) > 0;
            clause.push_back(alive ? -sat_var : sat_var);
        }
    }

    return clause;
}
