#pragma once
#include <set>
#include <string>

enum class SolverStatus {
    SAT,
    UNSAT,
    ERROR
};

struct SolverResult {
    // Default to UNSAT so an unconfigured/uninitialised SolverResult doesn't
    // get accidentally interpreted as a real SAT model with an empty solution.
    SolverStatus status = SolverStatus::UNSAT;
    std::set<int> solution;  // set of true literals (positive = true, negative = false)
    std::string error_message;
    std::set<int> failed_assumptions;
    std::set<int> core_variables;
};
