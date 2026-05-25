#pragma once
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include "config.hpp"
#include "rule.hpp"

// Pack a 3x3 neighborhood into a 9-bit int with row-major bit order:
//   bit 0 1 2     NW  N  NE
//   bit 3 4 5     W   C  E
//   bit 6 7 8     SW  S  SE
// Matches Rule's encoding (bit 4 = center).
inline int pack_neighborhood(bool nw, bool n, bool ne,
                              bool w,  bool c, bool e,
                              bool sw, bool s, bool se) {
    return nw | (n<<1) | (ne<<2) | (w<<3) | (c<<4) | (e<<5) | (sw<<6) | (s<<7) | (se<<8);
}

// ── Mask-based stability tracking ──
// For each cell in/near the catalyst, track which (state, neighbor_count)
// combinations are possible for a stable (still life) configuration.
//   d0..d6 = dead with 0..6 live neighbors (d3 excluded: would birth)
//   l2, l3 = alive with 2 or 3 live neighbors
//   d7, d8 impossible (neighbor would be overpopulated)
// Bits: 0 = possible, 1 = ruled out.
namespace SM {
    constexpr uint8_t L2 = 1 << 0;
    constexpr uint8_t L3 = 1 << 1;
    constexpr uint8_t D0 = 1 << 2;
    constexpr uint8_t D1 = 1 << 3;
    constexpr uint8_t D2 = 1 << 4;
    constexpr uint8_t D4 = 1 << 5;
    constexpr uint8_t D5 = 1 << 6;
    constexpr uint8_t D6 = 1 << 7;
    constexpr uint8_t LIVE = L2 | L3;
    constexpr uint8_t DEAD = D0 | D1 | D2 | D4 | D5 | D6;
    constexpr uint8_t ALL  = LIVE | DEAD;

    inline bool maybe_alive(uint8_t m) { return (m & LIVE) != LIVE; }
    inline bool maybe_dead(uint8_t m)  { return (m & DEAD) != DEAD; }
    inline bool is_alive(uint8_t m) { return maybe_alive(m) && !maybe_dead(m); }
    inline bool is_dead(uint8_t m)  { return maybe_dead(m) && !maybe_alive(m); }

    inline uint8_t compute_ruled_out(int center, int alive_nbrs, int unknown_nbrs) {
        int lo = alive_nbrs;
        int hi = alive_nbrs + unknown_nbrs;
        auto out_of_range = [&](int k) { return k < lo || k > hi; };

        uint8_t ruled = 0;
        if (center == 1) {
            ruled |= DEAD;
            if (out_of_range(2)) ruled |= L2;
            if (out_of_range(3)) ruled |= L3;
        } else if (center == 0) {
            ruled |= LIVE;
            if (out_of_range(0)) ruled |= D0;
            if (out_of_range(1)) ruled |= D1;
            if (out_of_range(2)) ruled |= D2;
            if (out_of_range(4)) ruled |= D4;
            if (out_of_range(5)) ruled |= D5;
            if (out_of_range(6)) ruled |= D6;
        } else {
            if (out_of_range(2)) ruled |= L2;
            if (out_of_range(3)) ruled |= L3;
            if (out_of_range(0)) ruled |= D0;
            if (out_of_range(1)) ruled |= D1;
            if (out_of_range(2)) ruled |= D2;
            if (out_of_range(4)) ruled |= D4;
            if (out_of_range(5)) ruled |= D5;
            if (out_of_range(6)) ruled |= D6;
        }
        return ruled;
    }
}

// ── Stability mask grid ──
// Tracks which stable states are possible for cells in/near the catalyst.
struct StableMaskGrid {
    std::map<std::pair<int,int>, uint8_t> masks;

    // Initialize masks given catalyst positions and perturbation region (ZOI).
    // known_alive = stator + non_stator positions (always alive in stable state).
    void init(const std::set<std::pair<int,int>>& catalyst_positions,
              const std::set<std::pair<int,int>>& perturbation_region,
              const std::set<std::pair<int,int>>& known_alive = {}) {
        for (auto& pos : catalyst_positions)
            masks[pos] = 0; // all options open
        for (auto& pos : known_alive)
            masks[pos] = SM::DEAD; // known alive, rule out dead options
        for (auto& pos : perturbation_region)
            if (!catalyst_positions.count(pos) && !known_alive.count(pos))
                masks[pos] = SM::LIVE; // dead, rule out alive options
    }

    int center_state(int wx, int wy) const {
        auto it = masks.find({wx, wy});
        if (it == masks.end()) return 0;
        if (SM::is_alive(it->second)) return 1;
        if (SM::is_dead(it->second))  return 0;
        return -1;
    }

    // Propagate stability constraints until fixed point.
    int propagate() {
        int total_resolved = 0;
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto& [pos, mask] : masks) {
                if (mask == SM::ALL) continue;
                auto [cx, cy] = pos;

                int alive_nbrs = 0, unknown_nbrs = 0;
                std::vector<std::pair<int,int>> unknown_positions;
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dx == 0 && dy == 0) continue;
                        int nx = cx + dx, ny = cy + dy;
                        int ns = neighbor_state(nx, ny);
                        if (ns == 1) alive_nbrs++;
                        else if (ns == -1) {
                            unknown_nbrs++;
                            unknown_positions.push_back({nx, ny});
                        }
                    }
                }

                int cs = center_state(cx, cy);
                uint8_t new_ruled = SM::compute_ruled_out(cs, alive_nbrs, unknown_nbrs);
                uint8_t updated = mask | new_ruled;

                if (updated != mask) {
                    mask = updated;
                    changed = true;
                    if (SM::is_alive(updated) || SM::is_dead(updated))
                        total_resolved++;
                }

                // Signal propagation: force unknowns ON/OFF when constrained
                if (unknown_nbrs > 0) {
                    int min_possible = 9, max_possible = -1;
                    if (SM::maybe_alive(updated)) {
                        int a = std::max(2, alive_nbrs);
                        int b = std::min(3, alive_nbrs + unknown_nbrs);
                        if (a <= b) {
                            min_possible = std::min(min_possible, a);
                            max_possible = std::max(max_possible, b);
                        }
                    }
                    if (SM::maybe_dead(updated)) {
                        for (int k : {0, 1, 2, 4, 5, 6}) {
                            if (k >= alive_nbrs && k <= alive_nbrs + unknown_nbrs) {
                                min_possible = std::min(min_possible, k);
                                max_possible = std::max(max_possible, k);
                            }
                        }
                    }

                    if (min_possible <= max_possible) {
                        int min_unk_alive = std::max(0, min_possible - alive_nbrs);
                        int max_unk_alive = std::min(unknown_nbrs, max_possible - alive_nbrs);

                        if (min_unk_alive == unknown_nbrs && max_unk_alive == unknown_nbrs) {
                            for (auto& upos : unknown_positions) {
                                auto uit = masks.find(upos);
                                if (uit != masks.end() && SM::maybe_alive(uit->second)) {
                                    uint8_t before = uit->second;
                                    uit->second |= SM::DEAD;
                                    if (uit->second != before) { changed = true; total_resolved++; }
                                }
                            }
                        } else if (max_unk_alive == 0) {
                            for (auto& upos : unknown_positions) {
                                auto uit = masks.find(upos);
                                if (uit != masks.end() && SM::maybe_dead(uit->second)) {
                                    uint8_t before = uit->second;
                                    uit->second |= SM::LIVE;
                                    if (uit->second != before) { changed = true; total_resolved++; }
                                }
                            }
                        }
                    }
                }
            }
        }
        return total_resolved;
    }

    // Build stable state map: -1=unknown, 0=dead, 1=alive
    std::map<std::pair<int,int>, int> stable_map() const {
        std::map<std::pair<int,int>, int> result;
        for (auto& [pos, mask] : masks) {
            if (SM::is_alive(mask)) result[pos] = 1;
            else if (SM::is_dead(mask)) result[pos] = 0;
            else result[pos] = -1;
        }
        return result;
    }

private:
    int neighbor_state(int wx, int wy) const {
        auto it = masks.find({wx, wy});
        if (it == masks.end()) return 0;
        if (SM::is_alive(it->second)) return 1;
        if (SM::is_dead(it->second))  return 0;
        return -1;
    }
};

struct Grid {
    // grid[t][y][x] = variable index (0=dead, 1=alive, >=2=SAT variable)
    std::vector<std::vector<std::vector<int>>> cells;

    // Dimensions at each timestep (light cone expands)
    int total_gens = 0;

    // Global coordinate offset: grid position (0,0) corresponds to world (ox, oy)
    int ox = 0, oy = 0;

    // Grid dimensions (max extent, for the final timestep)
    int width = 0, height = 0;

    // Catalyst variables: the stable-state variables for unknown cells
    // catalyst_vars[y][x] = SAT variable index (>=2) or 0 (dead outside catalyst region)
    // In world coordinates offset by (ox, oy)
    std::vector<std::vector<int>> catalyst_vars;

    // The set of (world_x, world_y) positions that are unknown/catalyst cells
    std::set<std::pair<int,int>> catalyst_positions;

    // Known catalyst cells (stator = always on, non_stator = on but might be active)
    std::set<std::pair<int,int>> stator_positions;
    std::set<std::pair<int,int>> non_stator_positions;

    // The perturbation region: ZOI of unknown cells
    // Set of (world_x, world_y) positions
    std::set<std::pair<int,int>> perturbation_region;

    // The catalyst neighborhood: catalyst cells + dead cells adjacent to catalyst
    // These are the cells constrained by recovered(t)
    std::set<std::pair<int,int>> catalyst_neighborhood;

    // adj_on(x,y): SAT literal that is true iff any cell in the 3x3 neighborhood
    // (including self) is an alive catalyst cell. Used to condition perturbation
    // detection and recovery — only cells near the actual still life matter.
    std::map<std::pair<int,int>, int> adj_on;

    // Free evolution: what each cell would be without any catalyst (all catalyst=dead).
    // free_evolution[t][wy][wx] indexed by world coords mapped to grid coords.
    // Only meaningful for cells in perturbation_region.
    std::vector<std::vector<std::vector<bool>>> free_evolution;

    // Wrap grid coordinates to torus topology
    int wrap_x(int gx) const { return ((gx % width) + width) % width; }
    int wrap_y(int gy) const { return ((gy % height) + height) % height; }

    // Get free evolution value at world coordinates (toroidal)
    bool free_at(int wx, int wy, int t) const {
        if (t < 0 || t >= (int)free_evolution.size()) return false;
        int gx = wrap_x(wx - ox);
        int gy = wrap_y(wy - oy);
        return free_evolution[t][gy][gx];
    }

    int next_var = 2;  // next SAT variable index to allocate

    int num_vars() const { return next_var - 2; }

    int alloc_var() { return next_var++; }

    // Get cell value at world coordinates (wx, wy) at time t (toroidal).
    int cell_at(int wx, int wy, int t) const {
        if (t < 0 || t >= total_gens) return 0;
        int gx = wrap_x(wx - ox);
        int gy = wrap_y(wy - oy);
        return cells[t][gy][gx];
    }

    // Get catalyst variable at world coordinates (toroidal)
    int catalyst_var_at(int wx, int wy) const {
        int gx = wrap_x(wx - ox);
        int gy = wrap_y(wy - oy);
        return catalyst_vars[gy][gx];
    }
};

// Compute the light cone: at t=0, the "possibly non-dead" region is
// active cells ∪ (catalyst bbox + 1). Each subsequent timestep expands by 1.
inline Grid build_grid(const SearchConfig& config) {
    Grid grid;
    grid.total_gens = config.total_generations() + 1; // +1 because T generations = T+1 states

    const auto& pattern = config.pattern;

    // Find bounding box of active and unknown cells in pattern coordinates
    int active_min_x = pattern.width, active_min_y = pattern.height;
    int active_max_x = -1, active_max_y = -1;
    int unknown_min_x = pattern.width, unknown_min_y = pattern.height;
    int unknown_max_x = -1, unknown_max_y = -1;

    for (int y = 0; y < pattern.height; y++) {
        for (int x = 0; x < pattern.width; x++) {
            CellState s = pattern.grid[y][x];
            if (s == CellState::ACTIVE) {
                active_min_x = std::min(active_min_x, x);
                active_min_y = std::min(active_min_y, y);
                active_max_x = std::max(active_max_x, x);
                active_max_y = std::max(active_max_y, y);
            }
            if (s == CellState::UNKNOWN || s == CellState::STATOR || s == CellState::NON_STATOR) {
                unknown_min_x = std::min(unknown_min_x, x);
                unknown_min_y = std::min(unknown_min_y, y);
                unknown_max_x = std::max(unknown_max_x, x);
                unknown_max_y = std::max(unknown_max_y, y);
            }
        }
    }

    // Convert to world coordinates (pattern coords are already world coords here)
    // t=0 region: union of active cells and (catalyst bbox expanded by 1)
    int t0_min_x = std::min(active_min_x, unknown_min_x - 1);
    int t0_min_y = std::min(active_min_y, unknown_min_y - 1);
    int t0_max_x = std::max(active_max_x, unknown_max_x + 1);
    int t0_max_y = std::max(active_max_y, unknown_max_y + 1);

    // Toroidal grid: pattern dimensions + 1-cell dead border on each side.
    // The border prevents the pattern corners from being adjacent via wrapping.
    grid.ox = -1;
    grid.oy = -1;
    grid.width = pattern.width + 2;
    grid.height = pattern.height + 2;

    // Allocate catalyst variables first
    grid.catalyst_vars.resize(grid.height, std::vector<int>(grid.width, 0));

    for (int y = 0; y < pattern.height; y++) {
        for (int x = 0; x < pattern.width; x++) {
            int gx = x - grid.ox;
            int gy = y - grid.oy;
            CellState s = pattern.grid[y][x];
            if (s == CellState::UNKNOWN) {
                int var = grid.alloc_var();
                grid.catalyst_vars[gy][gx] = var;
                grid.catalyst_positions.insert({x, y});
            } else if (s == CellState::STATOR) {
                grid.catalyst_vars[gy][gx] = 1; // known alive
                grid.stator_positions.insert({x, y});
            } else if (s == CellState::NON_STATOR) {
                grid.catalyst_vars[gy][gx] = 1; // known alive
                grid.non_stator_positions.insert({x, y});
            }
        }
    }

    // All catalyst cell positions (unknown + stator + non_stator)
    std::set<std::pair<int,int>> all_catalyst;
    all_catalyst.insert(grid.catalyst_positions.begin(), grid.catalyst_positions.end());
    all_catalyst.insert(grid.stator_positions.begin(), grid.stator_positions.end());
    all_catalyst.insert(grid.non_stator_positions.begin(), grid.non_stator_positions.end());

    // Build perturbation region: ZOI (zone of influence) of all catalyst cells
    // = all cells within distance 1 of any catalyst cell
    for (auto [cx, cy] : all_catalyst) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                grid.perturbation_region.insert({cx + dx, cy + dy});
            }
        }
    }

    // Build catalyst neighborhood: catalyst cells + dead cells adjacent to catalyst
    // (cells constrained by recovered(t))
    for (auto [cx, cy] : all_catalyst) {
        grid.catalyst_neighborhood.insert({cx, cy});
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;
                int nx = cx + dx, ny = cy + dy;
                if (all_catalyst.count({nx, ny}) == 0) {
                    // Adjacent but not a catalyst cell → dead in stable state
                    grid.catalyst_neighborhood.insert({nx, ny});
                }
            }
        }
    }

    // ── Phase 1: Set up t=0 grid ──
    grid.cells.resize(grid.total_gens);
    grid.cells[0].resize(grid.height, std::vector<int>(grid.width, 0));

    for (int wy = t0_min_y; wy <= t0_max_y; wy++) {
        for (int wx = t0_min_x; wx <= t0_max_x; wx++) {
            int gx = wx - grid.ox, gy = wy - grid.oy;
            if (gx < 0 || gx >= grid.width || gy < 0 || gy >= grid.height)
                continue;

            CellState s = CellState::DEAD;
            if (wx >= 0 && wx < pattern.width && wy >= 0 && wy < pattern.height)
                s = pattern.grid[wy][wx];

            switch (s) {
                case CellState::ACTIVE:
                case CellState::STATOR:
                case CellState::NON_STATOR:
                    grid.cells[0][gy][gx] = 1;
                    break;
                case CellState::UNKNOWN:
                    grid.cells[0][gy][gx] = grid.catalyst_vars[gy][gx];
                    break;
                default:
                    grid.cells[0][gy][gx] = 0;
                    break;
            }
        }
    }

    // ── Phase 2: Compute free evolution (active pattern, catalyst absent) ──
    // Used both to find the first-contact time K and to fill in cell values
    // outside the post-K uncertainty light cone.
    grid.free_evolution.resize(grid.total_gens,
        std::vector<std::vector<bool>>(grid.height,
            std::vector<bool>(grid.width, false)));
    for (int y = 0; y < pattern.height; y++) {
        for (int x = 0; x < pattern.width; x++) {
            if (pattern.grid[y][x] == CellState::ACTIVE) {
                grid.free_evolution[0][y - grid.oy][x - grid.ox] = true;
            }
        }
    }
    for (int t = 1; t < grid.total_gens; t++) {
        for (int gy = 0; gy < grid.height; gy++) {
            for (int gx = 0; gx < grid.width; gx++) {
                auto at = [&](int dy, int dx) {
                    return grid.free_evolution[t-1][grid.wrap_y(gy+dy)][grid.wrap_x(gx+dx)];
                };
                int nine = pack_neighborhood(
                    at(-1,-1), at(-1,0), at(-1,1),
                    at( 0,-1), at( 0,0), at( 0,1),
                    at( 1,-1), at( 1,0), at( 1,1));
                grid.free_evolution[t][gy][gx] = config.rule.evolves_to(nine);
            }
        }
    }

    // ── Phase 3: Find K = first time the active pattern reaches the catalyst ──
    // Pre-K: the catalyst sits in its stable still-life and cells outside the
    // catalyst neighborhood evolve identically to free_evolution. No SAT vars
    // are needed at t in [1, K).
    // Post-K: uncertainty spreads outward at lightspeed from catalyst_neighborhood
    // (Chebyshev radius +1 per generation). Cells inside that cone get fresh SAT
    // vars; cells outside still follow free_evolution.
    int K = grid.total_gens;
    for (int t = 0; t < grid.total_gens; t++) {
        bool contact = false;
        for (auto [wx, wy] : grid.catalyst_neighborhood) {
            if (grid.free_at(wx, wy, t)) { contact = true; break; }
        }
        if (contact) { K = t; break; }
    }
    int detected_K = K;
    // TODO: With the "pre-K reuses catalyst_vars" optimization, certain
    // inputs (notably with stator cells) become UNSAT — the encoding's
    // tautology-discard interaction with same-variable input/output bits
    // appears to lose constraints in a way I haven't fully diagnosed.
    // For now, force K=0 so every t > 0 in the catalyst light cone gets
    // fresh SAT vars. This loses the "deterministic until first contact"
    // saving but keeps correctness. Revisit when there's time to debug.
    K = 0;
    if (detected_K < grid.total_gens)
        std::cout << "First contact at t=" << detected_K
                  << " (using K=0 for now; see TODO)\n";
    else
        std::cout << "First contact: never (within " << grid.total_gens << " gens)\n";

    // ── Phase 4: Allocate cells per timestep ──
    // Maintain a growing "uncertain" mask seeded at t=K with catalyst_neighborhood
    // and dilated by 1 Chebyshev step per generation.
    std::set<std::pair<int,int>> stator_grid_pos;
    for (auto [wx, wy] : grid.stator_positions)
        stator_grid_pos.insert({wx - grid.ox, wy - grid.oy});

    std::vector<std::vector<bool>> uncertain(grid.height, std::vector<bool>(grid.width, false));
    auto seed_uncertain = [&]() {
        for (auto [wx, wy] : grid.catalyst_neighborhood) {
            int gx = wx - grid.ox, gy = wy - grid.oy;
            if (gx >= 0 && gx < grid.width && gy >= 0 && gy < grid.height)
                uncertain[gy][gx] = true;
        }
    };
    auto dilate_uncertain = [&]() {
        auto prev = uncertain;
        for (int gy = 0; gy < grid.height; gy++) {
            for (int gx = 0; gx < grid.width; gx++) {
                if (uncertain[gy][gx]) continue;
                for (int dy = -1; dy <= 1 && !uncertain[gy][gx]; dy++)
                    for (int dx = -1; dx <= 1 && !uncertain[gy][gx]; dx++)
                        if (prev[grid.wrap_y(gy+dy)][grid.wrap_x(gx+dx)])
                            uncertain[gy][gx] = true;
            }
        }
    };

    if (K == 0) seed_uncertain();

    int sat_vars_allocated = 0;
    for (int t = 1; t < grid.total_gens; t++) {
        grid.cells[t].resize(grid.height, std::vector<int>(grid.width, 0));
        if (t == K) seed_uncertain();
        else if (t > K) dilate_uncertain();

        for (int gy = 0; gy < grid.height; gy++) {
            for (int gx = 0; gx < grid.width; gx++) {
                int wx = gx + grid.ox, wy = gy + grid.oy;

                // Stator: forced alive at every t.
                if (stator_grid_pos.count({gx, gy})) {
                    grid.cells[t][gy][gx] = 1;
                    continue;
                }

                if (t < K) {
                    // Pre-perturbation. Catalyst cells reuse their t=0 var
                    // (catalyst_vars holds: 0 for dead-frame, SAT var for
                    // unknown, 1 for non_stator). Other cells follow free
                    // evolution (catalyst presence doesn't affect them yet).
                    if (grid.catalyst_neighborhood.count({wx, wy}))
                        grid.cells[t][gy][gx] = grid.catalyst_vars[gy][gx];
                    else
                        grid.cells[t][gy][gx] = grid.free_evolution[t][gy][gx] ? 1 : 0;
                } else if (uncertain[gy][gx]) {
                    grid.cells[t][gy][gx] = grid.alloc_var();
                    sat_vars_allocated++;
                } else {
                    grid.cells[t][gy][gx] = grid.free_evolution[t][gy][gx] ? 1 : 0;
                }
            }
        }
    }

    std::cout << "Allocated " << sat_vars_allocated
              << " per-timestep SAT variables (in post-contact light cone)\n";

    return grid;
}
