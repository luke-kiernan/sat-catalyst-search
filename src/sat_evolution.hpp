#pragma once
#include <vector>
#include <array>
#include <cassert>
#include "clause_types.hpp"

// Truth table for CGOL: 10-bit input (9 neighborhood + 1 output).
// Bit layout: bits 0-8 = 3x3 neighborhood in row-major order, bit 9 = next-gen output.
// table[x + 512*r] = 1 iff input x produces output r under CGOL rules.
static const std::vector<int> evolution_table = []() -> std::vector<int> {
    std::vector<int> table(1024);
    for (int x = 0; x < 512; x++) {
        int r = abs(__builtin_popcount(x & 495) * 2 + (x >> 4) % 2 - 6) <= 1;
        table[x + r * 512] = 1;
    }
    return table;
}();

// Prime implicants for NOT-CGOL over 10 variables.
// Each (care, force) pair: conjunction of all resulting clauses enforces CGOL.
static const std::vector<std::pair<int, int>> evolutionPrimeImplicants = []() {
    std::vector<std::pair<int, int>> ans;
    for (int care = 1; care < 1024; care++) {
        for (int force = care;; force = (force - 1) & care) {
            for (int x = care; x < 1024; x = (x + 1) | care)
                if (evolution_table[x ^ force])
                    goto no;
            for (auto [a, b] : ans)
                if ((a & care) == a && (a & force) == b)
                    goto no;
            ans.push_back({care, force});
        no:;
            if (force == 0)
                break;
        }
    }
    // Verify
    for (int i = 0; i < 1024; i++) {
        bool s = 1;
        for (auto [a, b] : ans)
            s = s && (a & (~(i ^ b)));
        assert(s == (bool)evolution_table[i]);
    }
    return ans;
}();

// Generate evolution clauses for a single transition.
// ten_cells[0..8] = neighborhood at time t (row-major), ten_cells[9] = output at t+1.
// Values: 0=dead, 1=alive, >=2 = SAT variable index.
inline PrimeClauseList generate_evolution_clauses(const std::array<int, 10>& ten_cells) {
    PrimeClauseList clauses;
    ClauseBuilder<10> builder;

    for (const auto& [care, force] : evolutionPrimeImplicants) {
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
