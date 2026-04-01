#pragma once
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <csignal>
#include <cstdlib>
#include "cadical.hpp"
#include "clause_types.hpp"
#include "solver_result.hpp"

// Global pointer for signal handler to reach the active solver
static CaDiCaL::Solver* g_active_solver = nullptr;
static struct sigaction g_old_sigint_action;

static void cadical_sigint_handler(int sig) {
    if (g_active_solver) {
        g_active_solver->terminate();
    } else {
        sigaction(SIGINT, &g_old_sigint_action, nullptr);
        raise(sig);
    }
}

class CadicalSolver {
public:
    CadicalSolver() : solver_(new CaDiCaL::Solver()) {
        solver_->set("quiet", 1);
        solver_->set("factor", 0);
    }

    ~CadicalSolver() {
        if (g_active_solver == solver_) {
            g_active_solver = nullptr;
        }
        delete solver_;
    }

    CadicalSolver(const CadicalSolver&) = delete;
    CadicalSolver& operator=(const CadicalSolver&) = delete;

    template<size_t N>
    void add_clause(const Clause<N>& clause) {
        for (int lit : clause) {
            if (lit != 0) solver_->add(lit);
        }
        solver_->add(0);
    }

    void add_clause(const std::vector<int>& clause) {
        for (int lit : clause) {
            solver_->add(lit);
        }
        solver_->add(0);
    }

    template<size_t N>
    void add_clauses(const ClauseList<N>& clauses) {
        for (const auto& clause : clauses) {
            add_clause(clause);
        }
    }

    void add_clauses(const BigClauseList& clauses) {
        for (const auto& clause : clauses) {
            add_clause(clause);
        }
    }

    void assume(int lit) {
        solver_->assume(lit);
    }

    SolverResult solve() {
        install_signal_handler();
        int status = solver_->solve();
        uninstall_signal_handler();
        return make_result(status);
    }

    int val(int var) const {
        return solver_->val(var);
    }

    CaDiCaL::Solver& raw() {
        return *solver_;
    }

private:
    CaDiCaL::Solver* solver_;

    void install_signal_handler() {
        g_active_solver = solver_;
        struct sigaction sa;
        sa.sa_handler = cadical_sigint_handler;
        sigemptyset(&sa.sa_mask);
        sa.sa_flags = 0;
        sigaction(SIGINT, &sa, &g_old_sigint_action);
    }

    void uninstall_signal_handler() {
        sigaction(SIGINT, &g_old_sigint_action, nullptr);
        g_active_solver = nullptr;
    }

    SolverResult make_result(int status) {
        SolverResult result;
        if (status == 10) {
            result.status = SolverStatus::SAT;
            int max_var = solver_->vars();
            for (int v = 1; v <= max_var; v++) {
                int val = solver_->val(v);
                result.solution.insert(val);
            }
        } else if (status == 20) {
            result.status = SolverStatus::UNSAT;
        } else {
            result.status = SolverStatus::ERROR;
            result.error_message = "Solver interrupted or reached limit";
        }
        return result;
    }
};
