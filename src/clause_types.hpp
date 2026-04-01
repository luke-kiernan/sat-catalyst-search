#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>

// Fixed-size clause: literals sorted, unused slots filled with sentinel 0
template<size_t N>
using Clause = std::array<int, N>;

template<size_t N>
using ClauseList = std::vector<Clause<N>>;

// Prime implicant encoding: max 9 literals (9-var stability) or 10 (10-var evolution)
using PrimeClause = Clause<10>;
using PrimeClauseList = ClauseList<10>;

// Stability clauses: max 9 literals
using StabilityClause = Clause<9>;
using StabilityClauseList = ClauseList<9>;

// For arbitrary-sized clauses
using BigClause = std::vector<int>;
using BigClauseList = std::vector<BigClause>;

// Helper for building clauses incrementally
// Detects tautologies (x and -x) and produces sorted Clause
template<size_t N>
class ClauseBuilder {
    Clause<N> lits{};
    int count = 0;
    bool tautology = false;

public:
    void clear() {
        for (int i = 0; i < count; i++) lits[i] = 0;
        count = 0;
        tautology = false;
    }

    bool add(int literal) {
        if (tautology) return true;
        for (int i = 0; i < count; i++) {
            if (lits[i] == -literal) {
                tautology = true;
                return true;
            }
        }
        if (count >= (int)N) {
            throw std::runtime_error("ClauseBuilder: clause exceeds max clause length");
        }
        lits[count++] = literal;
        return false;
    }

    bool is_tautology() const { return tautology; }
    bool empty() const { return count == 0; }

    Clause<N> get() const {
        Clause<N> c = lits;
        std::sort(c.begin(), c.begin() + count);
        return c;
    }
};

// Deduplicate a ClauseList in place
template<size_t N>
inline void deduplicate_clauses(ClauseList<N>& clauses) {
    std::sort(clauses.begin(), clauses.end());
    clauses.erase(std::unique(clauses.begin(), clauses.end()), clauses.end());
}
