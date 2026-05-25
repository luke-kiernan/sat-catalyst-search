#include <iostream>
#include <csignal>
#include "encoding.hpp"
#include "output.hpp"
#include "relevance.hpp"

// Find a minimal-population still life consistent with the relevant cells.
// Uses a second SAT solver with stability constraints + relevant cell fixings,
// binary searching on the number of non-relevant alive cells.
CatalystSolution minimize_catalyst(Grid& grid, const SolverResult& result,
                                    const std::set<std::pair<int,int>>& relevant) {
    // Collect non-relevant catalyst SAT variables
    std::vector<int> free_vars;
    for (auto [wx, wy] : grid.catalyst_positions) {
        if (relevant.count({wx, wy})) continue;
        int var = grid.catalyst_var_at(wx, wy);
        if (var >= 2)
            free_vars.push_back(var - 1);
    }

    int lo = 0, hi = (int)free_vars.size();
    SolverResult best;

    while (lo < hi) {
        int mid = (lo + hi) / 2;

        CadicalSolver attempt;
        encode_stability(attempt, grid);

        // Fix relevant cells
        for (auto [wx, wy] : relevant) {
            int var = grid.catalyst_var_at(wx, wy);
            if (var < 2) continue;
            int sat_var = var - 1;
            bool alive = result.solution.count(sat_var) > 0;
            attempt.add_clause(std::vector<int>{alive ? sat_var : -sat_var});
        }

        // At most mid non-relevant cells alive
        if (mid < (int)free_vars.size())
            encode_at_most_k(attempt, grid, free_vars, mid);

        SolverResult res = attempt.solve();
        if (res.status == SolverStatus::SAT) {
            best = res;
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }

    if (best.status != SolverStatus::SAT)
        return extract_catalyst(grid, result);

    return extract_catalyst(grid, best);
}

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
              << " non-stator\n";
    std::cout << "Perturbation region: " << grid.perturbation_region.size() << "\n\n";

    // Encode
    std::cout << "Encoding...\n";
    CadicalSolver solver;
    EncodingStats stats;
    TemporalVars tv = encode_all(solver, grid, config, stats);
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

        // Extract relevant cells and minimize catalyst
        solution_count++;
        auto relevant = find_relevant_cells(grid, result);
        CatalystSolution sol = minimize_catalyst(grid, result, relevant);
        std::cout << "--- Solution " << solution_count
                  << " (pop=" << sol.population
                  << ", pos=" << sol.min_x << "," << sol.min_y << ") ---\n";
        print_catalyst(sol);
        std::cout << catalyst_to_rle(sol) << "\n\n";
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
        std::string rle_file = write_summary_rle_file(all_solutions, argv[1]);
        if (!rle_file.empty())
            std::cout << "Results written to " << rle_file << "\n";
        std::string debug_file = write_debug_rle_file(all_debug_solutions, argv[1]);
        if (!debug_file.empty())
            std::cout << "Debug written to " << debug_file << "\n";
    }

    return 0;
}
