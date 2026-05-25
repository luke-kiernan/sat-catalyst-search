#include <iostream>
#include <csignal>
#include "encoding.hpp"
#include "output.hpp"
#include "relevance.hpp"

// NOTE: We report each solver solution as-is rather than minimizing its
// population. A previous "minimize" step fixed the relevant cells and searched
// for the smallest *still life* through them — but that is only sound when
// relevance provably determines the interaction. Our relevance heuristic
// ("cells adjacent to something that changes") makes no such guarantee, so
// stability-only minimization could strip alive cells that are structurally
// part of the catalyst, yielding a still life that no longer catalyzes.
// Reporting the solver's actual model avoids that entirely.

static volatile sig_atomic_t g_interrupted = 0;

static void sigint_handler(int) {
    g_interrupted = 1;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.toml>\n";
        return 1;
    }

    std::cout << "SAT Catalyst Search" << std::endl;
    std::cout << "Input: " << argv[1] << std::endl << std::endl;

    // Parse config
    std::cout << "Parsing config..." << std::endl;
    SearchConfig config = parse_config(argv[1]);
    std::cout << "Pattern: " << config.pattern.width << "x" << config.pattern.height
              << ", center=(" << config.center_x << "," << config.center_y << ")\n";
    std::cout << "Generations: " << config.total_generations() + 1 << " states" << std::endl << std::endl;

    // Build grid
    std::cout << "Building grid..." << std::endl;
    Grid grid = build_grid(config);
    std::cout << "Grid: " << grid.width << "x" << grid.height
              << ", origin=(" << grid.ox << "," << grid.oy << ")\n";
    std::cout << "Catalyst cells: " << grid.catalyst_positions.size()
              << " unknown, " << grid.stator_positions.size()
              << " stator, " << grid.non_stator_positions.size()
              << " non-stator, " << grid.init_off_positions.size()
              << " init-off\n";
    std::cout << "Perturbation region: " << grid.perturbation_region.size() << "\n\n";

    // Encode
    std::cout << "Encoding...\n";
    CadicalSolver solver;
    EncodingStats stats;
    RuleEncoding rule_enc = compute_rule_encoding(config.rule);
    TemporalVars tv = encode_all(solver, grid, config, stats, rule_enc);
    std::cout << "\n";

    // Install SIGINT handler for clean interruption
    struct sigaction sa;
    sa.sa_handler = sigint_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGINT, &sa, nullptr);

    // Incremental solve loop
    std::cout << "Solving...\n\n";
    int solution_count = 0;
    std::vector<FullSolution> all_solutions;       // minimized, for _results.rle
    std::vector<FullSolution> all_debug_solutions;  // original with relevance, for _debug.rle

    while (!g_interrupted) {
        SolverResult result = solver.solve();

        if (result.status == SolverStatus::UNSAT) {
            std::cout << "\nSearch complete (UNSAT). ";
            break;
        }

        if (result.status == SolverStatus::ERROR) {
            std::cout << "\nSolver interrupted.\n";
            break;
        }

        // Report the solver's catalyst as-is; relevance only drives the
        // blocking clause (search diversity) and the debug RLE.
        solution_count++;
        auto relevant = find_relevant_cells(grid, result);
        CatalystSolution sol = extract_catalyst(grid, result);
        std::cout << "--- Solution " << solution_count
                  << " (pop=" << sol.population
                  << ", pos=" << sol.min_x << "," << sol.min_y << ") ---\n";
        print_catalyst(sol);
        std::cout << catalyst_to_rle(sol, config.rule.canonical) << "\n\n";
        std::cout << "Relevant: " << relevant.size() << " / "
                  << grid.catalyst_positions.size() << " catalyst cells\n";

        std::vector<int> block;
        for (auto [wx, wy] : relevant) {
            int var = grid.catalyst_var_at(wx, wy);
            if (var < 2) continue;
            int sat_var = var - 1;
            bool alive = result.solution.count(sat_var) > 0;
            block.push_back(alive ? -sat_var : sat_var);
        }
        if (block.empty()) {
            block = blocking_clause(grid, result);
        }
        std::cout << "Blocking clause: " << block.size() << " literals\n\n";
        solver.add_clause(block);

        // Save minimized solution for main RLE, original for debug RLE
        all_solutions.push_back(extract_full_solution_from_catalyst(sol, config));
        all_debug_solutions.push_back(extract_full_solution(grid, result, config, &relevant));
    }

    std::cout << "Found " << solution_count << " solution(s).\n";

    // Write summary RLE (minimized) and debug RLE (original with relevance)
    if (!all_solutions.empty()) {
        std::string rle_file = write_summary_rle_file(all_solutions, argv[1], config.rule.canonical);
        if (!rle_file.empty())
            std::cout << "Results written to " << rle_file << "\n";
        std::string debug_file = write_debug_rle_file(all_debug_solutions, argv[1]);
        if (!debug_file.empty())
            std::cout << "Debug written to " << debug_file << "\n";
    }

    return 0;
}
