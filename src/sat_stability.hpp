#pragma once
#include <vector>
#include <array>
#include <cassert>
#include "clause_types.hpp"

// Truth table for "stable under CGOL": 9-bit input (center + 8 neighbors).
// Bit layout matches the 3x3 neighborhood scan order:
//   bit 0..8 = cells in row-major order (NW, N, NE, W, C, E, SW, S, SE)
//   bit 4 = center cell
// A configuration is stable iff:
//   center=1 and neighbor_count in {2,3}
//   center=0 and neighbor_count != 3
static const std::vector<int> stability_table = []() -> std::vector<int> {
    std::vector<int> table(512);
    for (int x = 0; x < 512; x++) {
        int center = (x >> 4) & 1;
        int neighbor_count = __builtin_popcount(x & ~(1 << 4));  // all bits except center
        bool stable;
        if (center == 1) {
            stable = (neighbor_count == 2 || neighbor_count == 3);
        } else {
            stable = (neighbor_count != 3);
        }
        table[x] = stable ? 1 : 0;
    }
    return table;
}();

// Prime implicants for NOT-stable-under-CGOL.
// Each (care, force) pair represents a minimal set of cell assignments that
// guarantees instability. Same algorithm as sat_logic.hpp but over 9 bits.
static const std::vector<std::pair<int, int>> stabilityPrimeImplicants = []() {
    std::vector<std::pair<int, int>> ans;
    for (int care = 1; care < 512; care++) {
        for (int force = care;; force = (force - 1) & care) {
            // Check: for all supersets of care with force applied, is the config unstable?
            for (int x = care; x < 512; x = (x + 1) | care)
                if (stability_table[x ^ force])  // if any matching config IS stable, not a valid implicant
                    goto no;
            // Minimality check: no existing implicant is a subset
            for (auto [a, b] : ans)
                if ((a & care) == a && (a & force) == b)
                    goto no;
            ans.push_back({care, force});
        no:;
            if (force == 0)
                break;
        }
    }
    // Verify: the conjunction of all prime implicant clauses exactly captures stability
    for (int i = 0; i < 512; i++) {
        bool all_clauses_satisfied = true;
        for (auto [a, b] : ans) {
            // Clause literal for care bit j: satisfied when cell_state == force_state
            // (cell matching force means it differs from the unstable AND-term requirement)
            bool clause_sat = false;
            for (int j = 0; j < 9; j++) {
                if (a & (1 << j)) {
                    bool cell_state = (i >> j) & 1;
                    bool force_state = (b >> j) & 1;
                    if (cell_state == force_state) {
                        clause_sat = true;
                        break;
                    }
                }
            }
            if (!clause_sat) {
                all_clauses_satisfied = false;
                break;
            }
        }
        assert(all_clauses_satisfied == (bool)stability_table[i]);
    }
    return ans;
}();

// Generate stability clauses for a single cell's 9-cell neighborhood.
// nine_cells[0..8] are variable indices (0=dead, 1=alive, >=2=SAT variable)
// in the same bit order as the truth table.
// Returns clauses that enforce stability (forbid all unstable configurations).
inline StabilityClauseList generate_stability_clauses(const std::array<int, 9>& nine_cells) {
    StabilityClauseList clauses;
    ClauseBuilder<9> builder;

    for (const auto& [care, force] : stabilityPrimeImplicants) {
        bool clause_satisfied = false;
        for (int bit = 0; bit < 9; bit++) {
            if (care & (1 << bit)) {
                int var_index = nine_cells[bit];
                if (var_index < 2) {  // known cell
                    bool force_state = (force >> bit) & 1;
                    bool cell_state = (var_index != 0);
                    if (cell_state == force_state) {
                        clause_satisfied = true;
                    }
                } else {
                    // unknown cell: literal = sign * SAT_var
                    int sign = (force & (1 << bit)) ? 1 : -1;
                    clause_satisfied = builder.add(sign * (var_index - 1));
                }
                if (clause_satisfied)
                    break;
            }
        }
        if (!clause_satisfied && !builder.empty())
            clauses.emplace_back(builder.get());
        builder.clear();
    }
    return clauses;
}
