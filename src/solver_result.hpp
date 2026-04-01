#pragma once
#include <set>
#include <string>

enum class SolverStatus {
    SAT,
    UNSAT,
    ERROR
};

struct SolverResult {
    SolverStatus status;
    std::set<int> solution;  // set of true literals (positive = true, negative = false)
    std::string error_message;
    std::set<int> failed_assumptions;
    std::set<int> core_variables;
};
