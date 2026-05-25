#pragma once
#include <vector>
#include <array>
#include <iostream>
#include <chrono>
#include "grid.hpp"
#include "sat_stability.hpp"
#include "sat_evolution.hpp"
#include "cadical_solver.hpp"

struct EncodingStats {
    int stability_clauses = 0;
    int evolution_clauses = 0;
    int temporal_clauses = 0;
    int perturbation_clauses = 0;
    int other_clauses = 0;

    void print() const {
        std::cout << "Clauses: stability=" << stability_clauses
                  << " evolution=" << evolution_clauses
                  << " temporal=" << temporal_clauses
                  << " perturbation=" << perturbation_clauses
                  << " other=" << other_clauses
                  << " total=" << (stability_clauses + evolution_clauses +
                                   temporal_clauses + perturbation_clauses + other_clauses)
                  << "\n";
    }
};

// Helper: add an implication a → b as clause (¬a ∨ b)
inline void add_implication(CadicalSolver& solver, int a, int b) {
    solver.add_clause(std::vector<int>{-a, b});
}

// Helper: add a → (b ∨ c) as clause (¬a ∨ b ∨ c)
inline void add_impl_or(CadicalSolver& solver, int a, int b, int c) {
    solver.add_clause(std::vector<int>{-a, b, c});
}

// Bundles the precomputed prime-implicant lists for a rule. Built once per
// run (compute_*_implicants is O(2^bits) and not free) and passed by const
// reference to the per-cell clause generators.
struct RuleEncoding {
    std::vector<std::pair<int,int>> stability_implicants;
    std::vector<std::pair<int,int>> evolution_implicants;
};

inline RuleEncoding compute_rule_encoding(const Rule& rule) {
    return {
        compute_stability_implicants(rule),
        compute_evolution_implicants(rule),
    };
}

// ── 1. Stability constraints ──────────────────────────────────────────────────
// For each catalyst cell (unknown + stator + non_stator), enforce that its
// 3x3 neighborhood is stable.
inline int encode_stability(CadicalSolver& solver, Grid& grid,
                            const RuleEncoding& enc) {
    int count = 0;

    // Collect all catalyst cells for stability constraints
    std::set<std::pair<int,int>> all_catalyst;
    all_catalyst.insert(grid.catalyst_positions.begin(), grid.catalyst_positions.end());
    all_catalyst.insert(grid.stator_positions.begin(), grid.stator_positions.end());
    all_catalyst.insert(grid.non_stator_positions.begin(), grid.non_stator_positions.end());

    for (auto [wx, wy] : all_catalyst) {
        std::array<int, 9> nine;
        int i = 0;
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                // Use catalyst variable if it's a catalyst cell, else 0 (dead)
                int nx = wx + dx, ny = wy + dy;
                nine[i++] = grid.catalyst_var_at(nx, ny);
            }
        }
        auto clauses = generate_stability_clauses(nine, enc.stability_implicants);
        for (const auto& c : clauses) {
            solver.add_clause(c);
        }
        count += clauses.size();
    }
    return count;
}

// ── 2. Evolution constraints ──────────────────────────────────────────────────
// For each cell in the light cone at each transition t→t+1.
inline int encode_evolution(CadicalSolver& solver, Grid& grid,
                            const RuleEncoding& enc) {
    int count = 0;
    for (int t = 0; t + 1 < grid.total_gens; t++) {
        for (int gy = 0; gy < grid.height; gy++) {
            for (int gx = 0; gx < grid.width; gx++) {
                // Output cell must be a variable (not boundary dead)
                int output = grid.cells[t + 1][gy][gx];
                if (output == 0) continue;  // dead boundary, no constraint needed

                int wx = gx + grid.ox;
                int wy = gy + grid.oy;

                std::array<int, 10> ten;
                int i = 0;
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        ten[i++] = grid.cell_at(wx + dx, wy + dy, t);
                    }
                }
                ten[9] = output;

                auto clauses = generate_evolution_clauses(ten, enc.evolution_implicants);
                for (const auto& c : clauses) {
                    solver.add_clause(c);
                }
                count += clauses.size();
            }
        }
    }
    return count;
}

// ── 3. Temporal tracking ──────────────────────────────────────────────────────
struct TemporalVars {
    // Per-generation variables
    std::vector<int> started;        // started(t): interaction has begun by gen t
    std::vector<int> recovered;      // recovered(t): catalyst matches stable state at gen t
    std::vector<int> perturbed;      // perturbed(t): at least one perturbation cell differs at gen t
    std::vector<int> must_recover;   // must_recover(t): recovery is required at gen t
};

inline TemporalVars encode_temporal(CadicalSolver& solver, Grid& grid,
                                     const SearchConfig& config, int& clause_count) {
    int T = grid.total_gens;
    TemporalVars tv;
    tv.started.resize(T);
    tv.recovered.resize(T);
    tv.perturbed.resize(T);
    tv.must_recover.resize(T);

    for (int t = 0; t < T; t++) {
        tv.started[t] = grid.alloc_var();
        tv.recovered[t] = grid.alloc_var();
        tv.perturbed[t] = grid.alloc_var();
        tv.must_recover[t] = grid.alloc_var();
    }

    clause_count = 0;

    // ── perturbed(t) definition ──
    // For catalyst cells: perturbed if cell(x,y,t) != catalyst_var(x,y)
    // For non-catalyst cells: perturbed if cell(x,y,t) != free_evolution(x,y,t)
    //   (detects birth inhibition: catalyst prevents a birth that free evolution would have)
    // Forward: diff(x,y,t) → perturbed(t)
    // Backward: perturbed(t) → ∨ diff(x,y,t)
    for (int t = 0; t < T; t++) {
        BigClause perturbed_definition;  // for the backward direction
        perturbed_definition.push_back(-tv.perturbed[t]);

        for (auto [wx, wy] : grid.perturbation_region) {
            // Stator cells never perturb (always alive by definition)
            // TODO(test): no direct coverage that perturbation is correctly
            // skipped for stators and constrained for non-stators (stable_val==1
            // path below). test_stator_encoding only checks end-to-end SAT.
            if (grid.stator_positions.count({wx, wy})) continue;

            int cell_t = grid.cell_at(wx, wy, t);
            int stable_val = grid.catalyst_var_at(wx, wy);

            if (stable_val >= 2) {
                // Unknown catalyst cell: perturbed if cell_t ≠ catalyst_var AND
                // cell is adjacent to an alive catalyst cell (adj_on).
                // Without the adj_on condition, dead catalyst cells would
                // create spurious perturbation when the glider passes through.
                if (cell_t >= 2 && cell_t != stable_val) {
                    auto adj_it = grid.adj_on.find({wx, wy});
                    if (adj_it != grid.adj_on.end()) {
                        int diff_var = grid.alloc_var();
                        int cell_lit = cell_t - 1;
                        int stable_lit = stable_val - 1;

                        // diff_var ↔ (cell_t ≠ stable_val) = XOR
                        solver.add_clause(std::vector<int>{-diff_var, cell_lit, stable_lit});
                        solver.add_clause(std::vector<int>{-diff_var, -cell_lit, -stable_lit});
                        solver.add_clause(std::vector<int>{-cell_lit, stable_lit, diff_var});
                        solver.add_clause(std::vector<int>{cell_lit, -stable_lit, diff_var});
                        clause_count += 4;

                        // Only counts as perturbation if adj_on:
                        // (diff_var ∧ adj_on) → perturbed(t)
                        int adj_diff = grid.alloc_var();
                        solver.add_clause(std::vector<int>{-adj_diff, diff_var});
                        solver.add_clause(std::vector<int>{-adj_diff, adj_it->second});
                        solver.add_clause(std::vector<int>{adj_diff, -diff_var, -adj_it->second});
                        clause_count += 3;

                        add_implication(solver, adj_diff, tv.perturbed[t]);
                        clause_count++;
                        perturbed_definition.push_back(adj_diff);
                    }
                }
            } else if (grid.non_stator_positions.count({wx, wy})) {
                // Non-stator catalyst cell: stable value = alive (1).
                // Perturbed if cell_t != 1 (cell turns off during interaction).
                // adj_on is trivially true (cell is itself a known-alive catalyst).
                // TODO: add a flag for whether non-stator cells count as perturbation
                // (for "interaction has started" and active cell counting).
                auto adj_it = grid.adj_on.find({wx, wy});
                if (adj_it == grid.adj_on.end()) continue;

                if (cell_t == 1) {
                    // Known alive = matches stable, no perturbation
                } else if (cell_t == 0) {
                    // Known dead ≠ stable alive → always perturbed (adj_on always true)
                    add_implication(solver, adj_it->second, tv.perturbed[t]);
                    clause_count++;
                    perturbed_definition.push_back(adj_it->second);
                } else {
                    // SAT variable: perturbed if cell is dead (not alive)
                    int diff_lit = -(cell_t - 1); // true when cell is dead
                    int adj_diff = grid.alloc_var();
                    solver.add_clause(std::vector<int>{-adj_diff, diff_lit});
                    solver.add_clause(std::vector<int>{-adj_diff, adj_it->second});
                    solver.add_clause(std::vector<int>{adj_diff, -diff_lit, -adj_it->second});
                    clause_count += 3;

                    add_implication(solver, adj_diff, tv.perturbed[t]);
                    clause_count++;
                    perturbed_definition.push_back(adj_diff);
                }
            } else {
                // Non-catalyst cell: perturbed if cell_t != free_evolution(x,y,t)
                // AND cell is adjacent to an alive catalyst cell (adj_on).
                auto adj_it = grid.adj_on.find({wx, wy});
                if (adj_it == grid.adj_on.end()) continue; // no catalyst neighbors

                bool free_val = grid.free_at(wx, wy, t);

                if (cell_t < 2) {
                    // Known cell: check against free evolution
                    bool actual = (cell_t == 1);
                    if (actual != free_val) {
                        // Perturbed if adj_on: (adj_on) → perturbed(t)
                        add_implication(solver, adj_it->second, tv.perturbed[t]);
                        clause_count++;
                        perturbed_definition.push_back(adj_it->second);
                    }
                } else {
                    // SAT variable: diff literal depends on free_val
                    int diff_lit;
                    if (free_val) {
                        diff_lit = -(cell_t - 1);
                    } else {
                        diff_lit = (cell_t - 1);
                    }
                    // active = adj_on ∧ diff
                    int adj_diff = grid.alloc_var();
                    solver.add_clause(std::vector<int>{-adj_diff, diff_lit});
                    solver.add_clause(std::vector<int>{-adj_diff, adj_it->second});
                    solver.add_clause(std::vector<int>{adj_diff, -diff_lit, -adj_it->second});
                    clause_count += 3;

                    add_implication(solver, adj_diff, tv.perturbed[t]);
                    clause_count++;
                    perturbed_definition.push_back(adj_diff);
                }
            }
        }

        // Backward: perturbed(t) → ∨ diff_lits
        solver.add_clause(perturbed_definition);
        clause_count++;
    }

    // ── started(t) definition ──
    // perturbed(t) → started(t)
    // started(t) → started(t+1)  (monotonicity)
    // started(t) → perturbed(t) ∨ started(t-1)
    for (int t = 0; t < T; t++) {
        add_implication(solver, tv.perturbed[t], tv.started[t]);
        clause_count++;

        if (t + 1 < T) {
            add_implication(solver, tv.started[t], tv.started[t + 1]);
            clause_count++;
        }

        if (t == 0) {
            // started(0) → perturbed(0)
            add_implication(solver, tv.started[0], tv.perturbed[0]);
            clause_count++;
        } else {
            // started(t) → perturbed(t) ∨ started(t-1)
            add_impl_or(solver, tv.started[t], tv.perturbed[t], tv.started[t - 1]);
            clause_count++;
        }
    }

    // ── recovered(t) definition ──
    // recovered(t) → for each (x,y) in catalyst neighborhood:
    //   ¬adj_on(x,y) ∨ (cell(x,y,t) = stable(x,y))
    // Cells far from any alive catalyst cell don't constrain recovery.
    //
    // NOTE: This is only the forward implication. The backward direction
    // (all cells match → recovered(t)) is omitted because recovered(t) is
    // only ever forced true from above (must_recover / min_stable_interval),
    // never derived from cell state. Adding the backward clauses would be
    // correct and could help solver propagation; omitted for now but may
    // matter if recovered(t) is ever used where the solver needs to infer it.
    for (int t = 0; t < T; t++) {
        for (auto [wx, wy] : grid.catalyst_neighborhood) {
            // Stator cells: always alive, always match stable → skip
            if (grid.stator_positions.count({wx, wy})) continue;

            auto adj_it = grid.adj_on.find({wx, wy});
            if (adj_it == grid.adj_on.end()) continue;  // no catalyst neighbors
            int adj_lit = adj_it->second;

            int cell_t = grid.cell_at(wx, wy, t);
            int stable_val = grid.catalyst_var_at(wx, wy);

            if (stable_val >= 2) {
                // Unknown catalyst cell: recovered(t) → ¬adj_on ∨ (cell_t ↔ stable_val)
                if (cell_t >= 2 && cell_t != stable_val) {
                    int cell_lit = cell_t - 1;
                    int stable_lit = stable_val - 1;
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, -cell_lit, stable_lit});
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, cell_lit, -stable_lit});
                    clause_count += 2;
                } else if (cell_t == 0) {
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, -(stable_val - 1)});
                    clause_count++;
                } else if (cell_t == 1) {
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, stable_val - 1});
                    clause_count++;
                }
            } else if (grid.non_stator_positions.count({wx, wy})) {
                // Non-stator catalyst cell: stable = alive (1).
                // recovered(t) → ¬adj_on ∨ cell_t is alive
                // (adj_on is trivially true here since the cell itself is known-alive,
                // but we keep it in the clause and let unit propagation simplify)
                if (cell_t >= 2) {
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, cell_t - 1});
                    clause_count++;
                } else if (cell_t == 0) {
                    // Known dead ≠ stable alive: recovery impossible if adj_on
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit});
                    clause_count++;
                }
                // cell_t == 1: matches stable, no constraint needed
            } else {
                // Non-catalyst cell (stable_val=0, dead in stable state):
                // recovered(t) → ¬adj_on ∨ cell_t is dead
                if (cell_t >= 2) {
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit, -(cell_t - 1)});
                    clause_count++;
                } else if (cell_t == 1) {
                    // Known alive, must have ¬adj_on for recovery
                    solver.add_clause(std::vector<int>{-tv.recovered[t], -adj_lit});
                    clause_count++;
                }
            }
        }
    }

    // ── must_recover(t) ──
    // perturbed(t) → must_recover(t + max_active_window)
    int max_window = config.active_window_range[1];
    for (int t = 0; t < T; t++) {
        int recovery_t = t + max_window;
        if (recovery_t < T) {
            add_implication(solver, tv.perturbed[t], tv.must_recover[recovery_t]);
            clause_count++;
        }
    }
    // must_recover(t) → recovered(t)
    for (int t = 0; t < T; t++) {
        add_implication(solver, tv.must_recover[t], tv.recovered[t]);
        clause_count++;
    }

    // ── 4. Timing constraints ──

    // first-active-range [a, b]: ¬started(a-1) and started(b)
    int first_a = config.first_active_range[0];
    int first_b = config.first_active_range[1];
    if (first_a > 0 && first_a - 1 < T) {
        solver.add_clause(std::vector<int>{-tv.started[first_a - 1]});
        clause_count++;
    }
    if (first_b < T) {
        solver.add_clause(std::vector<int>{tv.started[first_b]});
        clause_count++;
    }

    // min-stable-interval: last min_stable_interval generations must be recovered
    int stable_start = T - config.min_stable_interval;
    for (int t = stable_start; t < T; t++) {
        if (t >= 0) {
            solver.add_clause(std::vector<int>{tv.recovered[t]});
            clause_count++;
        }
    }

    return tv;
}

// ── 5. Non-triviality ─────────────────────────────────────────────────────────
inline int encode_nontriviality(CadicalSolver& solver, const Grid& grid) {
    // If stator/non_stator cells exist, catalyst is already non-trivial
    if (!grid.stator_positions.empty() || !grid.non_stator_positions.empty())
        return 0;

    BigClause at_least_one;
    for (auto [wx, wy] : grid.catalyst_positions) {
        int var = grid.catalyst_var_at(wx, wy);
        if (var >= 2) {
            at_least_one.push_back(var - 1);
        }
    }
    solver.add_clause(at_least_one);
    return 1;
}

// ── 6. Perturbation constraints ───────────────────────────────────────────────

// Sequential counter encoding for at-most-k constraint.
// Given n literals, creates O(n*k) clauses and O(n*k) auxiliary variables.
inline int encode_at_most_k(CadicalSolver& solver, Grid& grid,
                             const std::vector<int>& lits, int k) {
    int n = lits.size();
    if (k >= n) return 0;  // trivially satisfied
    if (k < 0) return 0;

    int clause_count = 0;

    // Counter variables: counter[i][j] = "at least j+1 of lits[0..i] are true"
    // i in [0, n-1], j in [0, k]
    std::vector<std::vector<int>> counter(n, std::vector<int>(k + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= k; j++) {
            counter[i][j] = grid.alloc_var();
        }
    }

    // Base case: counter[0][0] ↔ lits[0]
    add_implication(solver, lits[0], counter[0][0]);
    add_implication(solver, counter[0][0], lits[0]);
    clause_count += 2;

    // counter[0][j] = false for j > 0
    for (int j = 1; j <= k; j++) {
        solver.add_clause(std::vector<int>{-(counter[0][j])});
        clause_count++;
    }

    // Inductive case
    for (int i = 1; i < n; i++) {
        // counter[i][0] ↔ lits[i] ∨ counter[i-1][0]
        add_implication(solver, lits[i], counter[i][0]);
        add_implication(solver, counter[i - 1][0], counter[i][0]);
        solver.add_clause(std::vector<int>{-counter[i][0], lits[i], counter[i - 1][0]});
        clause_count += 3;

        for (int j = 1; j <= k; j++) {
            // counter[i][j] ↔ (lits[i] ∧ counter[i-1][j-1]) ∨ counter[i-1][j]
            // Forward:
            solver.add_clause(std::vector<int>{-lits[i], -counter[i - 1][j - 1], counter[i][j]});
            add_implication(solver, counter[i - 1][j], counter[i][j]);
            clause_count += 2;

            // Backward:
            solver.add_clause(std::vector<int>{-counter[i][j], lits[i], counter[i - 1][j]});
            solver.add_clause(std::vector<int>{-counter[i][j], counter[i - 1][j - 1], counter[i - 1][j]});
            clause_count += 2;
        }
    }

    // At most k: ¬counter[n-1][k]
    solver.add_clause(std::vector<int>{-(counter[n - 1][k])});
    clause_count++;

    return clause_count;
}

// active(x,y,t): cell differs from its stable value, conditioned on adj_on.
// Only counts as active if the cell is adjacent to an alive catalyst cell.
// For non-catalyst cells: stable_val=0 (dead), so active = adj_on ∧ (cell alive)
// For unknown catalyst cells: active = adj_on ∧ (cell_t ≠ stable_val)
// For non-stator catalyst cells: active = adj_on ∧ (cell_t ≠ 1) = adj_on ∧ ¬cell_t
// Stator cells are skipped (never active).
inline int get_active_literal(Grid& grid, CadicalSolver& solver, int wx, int wy, int t,
                               int& clause_count) {
    // Stator cells are never active
    if (grid.stator_positions.count({wx, wy})) return 0;

    int cell_t = grid.cell_at(wx, wy, t);
    int stable_val = grid.catalyst_var_at(wx, wy);

    if (cell_t < 2) return 0;  // known cell, can't be "active" variable

    auto adj_it = grid.adj_on.find({wx, wy});
    if (adj_it == grid.adj_on.end()) return 0;  // no catalyst neighbors → never active
    int adj_lit = adj_it->second;

    int raw_diff;
    if (stable_val == 0) {
        // Dead in stable state → differs iff alive at t
        raw_diff = cell_t - 1;
    } else if (stable_val == 1) {
        // Non-stator: alive in stable state → differs iff dead at t
        raw_diff = -(cell_t - 1);
    } else if (cell_t == stable_val) {
        return 0;  // same variable, always equal
    } else {
        // Unknown catalyst var: differs iff cell_t ≠ stable_val (XOR)
        int diff_var = grid.alloc_var();
        int cell_lit = cell_t - 1;
        int stable_lit = stable_val - 1;
        solver.add_clause(std::vector<int>{-diff_var, cell_lit, stable_lit});
        solver.add_clause(std::vector<int>{-diff_var, -cell_lit, -stable_lit});
        solver.add_clause(std::vector<int>{-cell_lit, stable_lit, diff_var});
        solver.add_clause(std::vector<int>{cell_lit, -stable_lit, diff_var});
        clause_count += 4;
        raw_diff = diff_var;
    }

    // active = adj_on ∧ raw_diff
    int active_var = grid.alloc_var();
    solver.add_clause(std::vector<int>{-active_var, raw_diff});
    solver.add_clause(std::vector<int>{-active_var, adj_lit});
    solver.add_clause(std::vector<int>{active_var, -raw_diff, -adj_lit});
    clause_count += 3;

    return active_var;
}

inline int encode_perturbation(CadicalSolver& solver, Grid& grid,
                                const SearchConfig& config) {
    int clause_count = 0;
    int T = grid.total_gens;

    // max-active-cells: at-most-k active cells per generation
    if (config.max_active_cells >= 0) {
        for (int t = 0; t < T; t++) {
            std::vector<int> active_lits;
            for (auto [wx, wy] : grid.perturbation_region) {
                int lit = get_active_literal(grid, solver, wx, wy, t, clause_count);
                if (lit != 0) active_lits.push_back(lit);
            }
            if (!active_lits.empty()) {
                clause_count += encode_at_most_k(solver, grid, active_lits,
                                                  config.max_active_cells);
            }
        }
    }

    // max-ever-active-cells: at-most-k cells that are ever active
    if (config.max_ever_active_cells >= 0) {
        std::vector<int> ever_active_lits;
        for (auto [wx, wy] : grid.perturbation_region) {
            // ever_active(x,y) = ∨_t active(x,y,t)
            std::vector<int> per_time_lits;
            for (int t = 0; t < T; t++) {
                int lit = get_active_literal(grid, solver, wx, wy, t, clause_count);
                if (lit != 0) per_time_lits.push_back(lit);
            }
            if (!per_time_lits.empty()) {
                int ever_var = grid.alloc_var();
                // active(x,y,t) → ever_active(x,y)
                for (int lit : per_time_lits) {
                    add_implication(solver, lit, ever_var);
                    clause_count++;
                }
                // ever_active → ∨ active(x,y,t)
                BigClause backward;
                backward.push_back(-ever_var);
                for (int lit : per_time_lits) backward.push_back(lit);
                solver.add_clause(backward);
                clause_count++;

                ever_active_lits.push_back(ever_var);
            }
        }
        if (!ever_active_lits.empty()) {
            clause_count += encode_at_most_k(solver, grid, ever_active_lits,
                                              config.max_ever_active_cells);
        }
    }

    // ever-active-bounds [w, h]: pairwise exclusion for cells too far apart
    if (config.ever_active_bounds[0] >= 0 && config.ever_active_bounds[1] >= 0) {
        int bw = config.ever_active_bounds[0];
        int bh = config.ever_active_bounds[1];

        // Build ever_active variables for each perturbation cell
        // (reuse if already built above, but simpler to just create new ones)
        std::vector<std::pair<std::pair<int,int>, int>> cell_vars;  // ((wx,wy), ever_active_var)
        for (auto [wx, wy] : grid.perturbation_region) {
            std::vector<int> per_time_lits;
            for (int t = 0; t < T; t++) {
                int lit = get_active_literal(grid, solver, wx, wy, t, clause_count);
                if (lit != 0) per_time_lits.push_back(lit);
            }
            if (!per_time_lits.empty()) {
                int ever_var = grid.alloc_var();
                for (int lit : per_time_lits) {
                    add_implication(solver, lit, ever_var);
                    clause_count++;
                }
                BigClause backward;
                backward.push_back(-ever_var);
                for (int lit : per_time_lits) backward.push_back(lit);
                solver.add_clause(backward);
                clause_count++;

                cell_vars.push_back({{wx, wy}, ever_var});
            }
        }

        // Pairwise: if two cells are too far apart, they can't both be ever-active
        for (size_t i = 0; i < cell_vars.size(); i++) {
            for (size_t j = i + 1; j < cell_vars.size(); j++) {
                auto [pos_i, var_i] = cell_vars[i];
                auto [pos_j, var_j] = cell_vars[j];
                int dx = std::abs(pos_i.first - pos_j.first);
                int dy = std::abs(pos_i.second - pos_j.second);
                if (dx >= bw || dy >= bh) {
                    solver.add_clause(std::vector<int>{-var_i, -var_j});
                    clause_count++;
                }
            }
        }
    }

    return clause_count;
}

// ── adj_on computation ───────────────────────────────────────────────────────
// Populates grid.adj_on: for each cell in catalyst_neighborhood,
// adj_on = OR of catalyst vars in its 3x3 neighborhood.
// Stator/non_stator cells (catalyst_var == 1) make adj_on trivially true.
inline int encode_adj_on(CadicalSolver& solver, Grid& grid) {
    int clause_count = 0;

    // Allocate a single always-true variable for cells with known-alive neighbors
    int true_var = grid.alloc_var();
    solver.add_clause(std::vector<int>{true_var});
    clause_count++;

    for (auto [wx, wy] : grid.catalyst_neighborhood) {
        bool has_known_alive = false;
        std::vector<int> neighbor_lits;
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                int nvar = grid.catalyst_var_at(wx + dx, wy + dy);
                if (nvar == 1) {
                    has_known_alive = true; // stator or non_stator neighbor
                } else if (nvar >= 2) {
                    neighbor_lits.push_back(nvar - 1);
                }
            }
        }

        if (has_known_alive) {
            // Trivially true: at least one neighbor is always alive
            grid.adj_on[{wx, wy}] = true_var;
        } else if (neighbor_lits.empty()) {
            continue;
        } else if (neighbor_lits.size() == 1) {
            grid.adj_on[{wx, wy}] = neighbor_lits[0];
        } else {
            int adj_var = grid.alloc_var();
            for (int lit : neighbor_lits) {
                add_implication(solver, lit, adj_var);
                clause_count++;
            }
            BigClause backward;
            backward.push_back(-adj_var);
            for (int lit : neighbor_lits) backward.push_back(lit);
            solver.add_clause(backward);
            clause_count++;
            grid.adj_on[{wx, wy}] = adj_var;
        }
    }
    return clause_count;
}

// ── Main encoding entry point ─────────────────────────────────────────────────
inline TemporalVars encode_all(CadicalSolver& solver, Grid& grid,
                                const SearchConfig& config, EncodingStats& stats,
                                const RuleEncoding& rule_enc) {
    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "  Rule '" << config.rule.canonical << "': "
              << rule_enc.stability_implicants.size() << " stability + "
              << rule_enc.evolution_implicants.size() << " evolution implicants\n";

    encode_adj_on(solver, grid);

    stats.stability_clauses = encode_stability(solver, grid, rule_enc);
    std::cout << "  Stability: " << stats.stability_clauses << " clauses" << std::endl;

    stats.evolution_clauses = encode_evolution(solver, grid, rule_enc);
    std::cout << "  Evolution: " << stats.evolution_clauses << " clauses" << std::endl;

    int temporal_count = 0;
    TemporalVars tv = encode_temporal(solver, grid, config, temporal_count);
    stats.temporal_clauses = temporal_count;
    std::cout << "  Temporal: " << stats.temporal_clauses << " clauses" << std::endl;

    stats.perturbation_clauses = encode_perturbation(solver, grid, config);
    std::cout << "  Perturbation: " << stats.perturbation_clauses << " clauses" << std::endl;

    stats.other_clauses = encode_nontriviality(solver, grid);
    std::cout << "  Non-triviality: " << stats.other_clauses << " clauses" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Total encoding: " << ms << " ms, " << grid.num_vars() << " variables\n";
    stats.print();

    return tv;
}
