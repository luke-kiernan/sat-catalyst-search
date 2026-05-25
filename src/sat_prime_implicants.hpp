#pragma once
#include <vector>
#include <utility>
#include <cassert>

// Compute prime implicants of the COMPLEMENT of a boolean truth table.
//
// table has size 2^num_bits; table[i] = 1 means input i is in the positive set.
// Each returned (care, force) pair encodes a CNF clause that excludes one
// minimal subset of the negative set: literal for care bit j is satisfied
// when input bit j differs from force bit j.
//
// Identical algorithm for stability (num_bits=9) and evolution (num_bits=10);
// the two prior call sites differed only in truth table and bit width.
inline std::vector<std::pair<int,int>> compute_prime_implicants(
    const std::vector<int>& table, int num_bits)
{
    int sz = 1 << num_bits;
    assert((int)table.size() == sz);
    std::vector<std::pair<int,int>> ans;

    for (int care = 1; care < sz; care++) {
        for (int force = care;; force = (force - 1) & care) {
            // All inputs matching (care, force) must be in the negative set.
            for (int x = care; x < sz; x = (x + 1) | care)
                if (table[x ^ force]) goto skip;
            // Minimality: no existing implicant is a subset.
            for (auto [a, b] : ans)
                if ((a & care) == a && (a & force) == b) goto skip;
            ans.push_back({care, force});
        skip:;
            if (force == 0) break;
        }
    }

    // Verify: AND of all implicant clauses exactly equals the table.
    for (int i = 0; i < sz; i++) {
        bool s = true;
        for (auto [a, b] : ans)
            s = s && (a & (~(i ^ b)));
        assert(s == (bool)table[i]);
    }
    return ans;
}
