#pragma once
#include <vector>
#include <array>
#include "clause_types.hpp"
#include "rule.hpp"
#include "sat_prime_implicants.hpp"

// Build a 10-bit "is-correct-transition" truth table for the given rule.
// Bits 0..8 = 3x3 neighborhood (row-major, bit 4 = center at time t),
// bit 9 = output at time t+1.
// table[n + 512*r] = 1 iff input neighborhood n produces output r under the rule.
inline std::vector<int> build_evolution_table(const Rule& rule) {
    std::vector<int> table(1024);
    for (int n = 0; n < 512; n++) {
        int r = rule.evolves_to(n) ? 1 : 0;
        table[n + 512 * r] = 1;
    }
    return table;
}

inline std::vector<std::pair<int,int>> compute_evolution_implicants(const Rule& rule) {
    return compute_prime_implicants(build_evolution_table(rule), 10);
}

// Generate evolution clauses for a single transition under a precomputed
// implicant list. ten_cells[0..8] = neighborhood at time t (row-major),
// ten_cells[9] = output at time t+1.
inline PrimeClauseList generate_evolution_clauses(
    const std::array<int, 10>& ten_cells,
    const std::vector<std::pair<int,int>>& implicants)
{
    PrimeClauseList clauses;
    ClauseBuilder<10> builder;

    for (const auto& [care, force] : implicants) {
        bool clause_satisfied = false;
        for (int bit = 0; bit < 10; bit++) {
            if (care & (1 << bit)) {
                int var_index = ten_cells[bit];
                if (var_index < 2) {  // known cell
                    bool force_state = (force & (1 << bit)) != 0;
                    bool cell_state = (var_index != 0);
                    if (cell_state == force_state)
                        clause_satisfied = true;
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
