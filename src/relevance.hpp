#pragma once
#include <set>
#include "grid.hpp"
#include "solver_result.hpp"

// Determine which catalyst cells are "relevant" to the catalytic interaction.
//
// Rule-agnostic algorithm: resolve this solution's full evolution (every cell
// at every generation), then find every cell that *changes state* between two
// consecutive generations. Any UNKNOWN catalyst cell within the 3x3
// neighborhood of such a changing cell is marked relevant.
//
// Rationale: a cell only changes if it's part of the active pattern or has been
// perturbed by it, so the changing cells trace out exactly the interaction. A
// catalyst cell adjacent to that interaction is doing (or could be doing) the
// catalytic work; one that never neighbors a change is an untouched part of the
// still life and can be minimized away. This deliberately over-approximates
// (flags more than strictly necessary) in exchange for being simple and
// rule-independent — no Game-of-Life-specific neighbor thresholds.

inline std::set<std::pair<int,int>> find_relevant_cells(
    const Grid& grid, const SolverResult& result)
{
    int T = grid.total_gens;
    int H = grid.height, W = grid.width;

    // Resolve a grid cell to its concrete state under this solution.
    // grid.cells holds 0 (dead), 1 (alive), or a SAT variable (>=2).
    auto resolve = [&](int t, int gy, int gx) -> int {
        int v = grid.cells[t][gy][gx];
        if (v <= 1) return v;
        return result.solution.count(v - 1) > 0 ? 1 : 0;
    };

    std::set<std::pair<int,int>> relevant;

    for (int t = 0; t + 1 < T; t++) {
        for (int gy = 0; gy < H; gy++) {
            for (int gx = 0; gx < W; gx++) {
                if (resolve(t, gy, gx) == resolve(t + 1, gy, gx)) continue;
                // This cell changes state — flag the catalyst cells around it.
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        int ny = grid.wrap_y(gy + dy);
                        int nx = grid.wrap_x(gx + dx);
                        int wx = nx + grid.ox, wy = ny + grid.oy;
                        if (grid.catalyst_positions.count({wx, wy}))
                            relevant.insert({wx, wy});
                    }
                }
            }
        }
    }

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

    // TODO(test): the empty-relevant fallback is uncovered. Hard to trigger
    // organically — would need a non-empty catalyst where forward-sim resolves
    // every cell without peeking. Construct synthetically with a hand-built
    // Grid + SolverResult rather than from a solved instance.
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
