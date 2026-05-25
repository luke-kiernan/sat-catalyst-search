#pragma once
#include <vector>
#include <array>
#include "clause_types.hpp"
#include "rule.hpp"
#include "sat_prime_implicants.hpp"

// Build a 9-bit "is-stable" truth table for the given rule.
// Bit layout (matches 3x3 row-major scan): bit 4 = center, bits 0..8 = NW,N,NE,W,C,E,SW,S,SE.
// table[n] = 1 iff configuration n is a still life under the rule.
inline std::vector<int> build_stability_table(const Rule& rule) {
    std::vector<int> table(512);
    for (int n = 0; n < 512; n++)
        table[n] = rule.is_stable(n) ? 1 : 0;
    return table;
}

inline std::vector<std::pair<int,int>> compute_stability_implicants(const Rule& rule) {
    return compute_prime_implicants(build_stability_table(rule), 9);
}

// Generate stability clauses for a single cell's 9-cell neighborhood under
// a precomputed implicant list. nine_cells[0..8] are variable indices in
// row-major order (0=dead, 1=alive, >=2=SAT variable).
inline StabilityClauseList generate_stability_clauses(
    const std::array<int, 9>& nine_cells,
    const std::vector<std::pair<int,int>>& implicants)
{
    StabilityClauseList clauses;
    ClauseBuilder<9> builder;

    for (const auto& [care, force] : implicants) {
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
                    int sign = (force & (1 << bit)) ? 1 : -1;
                    clause_satisfied = builder.add(sign * (var_index - 1));
                }
                if (clause_satisfied) break;
            }
        }
        if (!clause_satisfied && !builder.empty())
            clauses.emplace_back(builder.get());
        builder.clear();
    }
    return clauses;
}
