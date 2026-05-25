#include <iostream>
#include <cassert>
#include <cstdlib>
#include "config.hpp"
#include "grid.hpp"
#include "encoding.hpp"
#include "relevance.hpp"
#include "cadical_solver.hpp"

// All three tests share one solved instance — encoding + solve dominates runtime.
int main() {
    std::cout << "=== Relevance / smart blocking tests ===\n";

    SearchConfig cfg = parse_config("inputs/glider_tight.toml");
    Grid grid = build_grid(cfg);

    CadicalSolver solver;
    EncodingStats stats;
    RuleEncoding rule_enc = compute_rule_encoding(cfg.rule);
    TemporalVars tv = encode_all(solver, grid, cfg, stats, rule_enc);
    (void)tv;

    SolverResult result = solver.solve();
    assert(result.status == SolverStatus::SAT);

    // 1. Smart blocking actually narrows the blocking set: the relevant
    //    set must be a strict subset of the full catalyst region.
    //    (Without this, the "smart" path could degenerate to naive blocking
    //    and the other assertions would still pass.)
    std::cout << "-- find_relevant_cells narrows the blocking set\n";
    auto relevant = find_relevant_cells(grid, result);
    std::cout << "   relevant=" << relevant.size()
              << " / catalyst=" << grid.catalyst_positions.size() << "\n";
    assert(!relevant.empty());
    assert(relevant.size() < grid.catalyst_positions.size());

    // 2. Clause polarity: each literal must be the negation of its variable's
    //    truth value in result.solution. (Trivially-easy line to invert.)
    // 3. Blocking correctness: every literal in the clause evaluates to FALSE
    //    under result.solution → adding the clause forbids exactly this
    //    assignment. This is the central correctness contract.
    std::cout << "-- smart_blocking_clause has correct polarity and blocks the solution\n";
    auto clause = smart_blocking_clause(grid, result);
    assert(!clause.empty());
    for (int lit : clause) {
        int sat_var = std::abs(lit);
        bool var_true = result.solution.count(sat_var) > 0;
        // lit is satisfied iff: (lit>0 && var_true) || (lit<0 && !var_true).
        // For the clause to block the solution, every literal must be unsatisfied.
        bool lit_satisfied = (lit > 0) ? var_true : !var_true;
        assert(!lit_satisfied);
        // Polarity check: positive lit ↔ variable was false in solution.
        assert((lit > 0) == !var_true);
    }
    std::cout << "   clause size=" << clause.size() << "\n";

    std::cout << "=== All relevance tests passed ===\n";
    return 0;
}
