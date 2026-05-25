#include <iostream>
#include <cassert>
#include "config.hpp"
#include "grid.hpp"
#include "encoding.hpp"
#include "cadical_solver.hpp"
#include "output.hpp"

// Integration test: encode + solve a stator-containing input end to end.
// Catches regressions in the encoding paths that special-case stator and
// non-stator cells (stability over union of catalyst sets, perturbation
// skipping stators, recovery treating stators/non-stators correctly,
// get_active_literal for known-alive stable values, encode_adj_on nvar==1).
static void run_case(const std::string& path,
                     int stator_x, int stator_y,
                     bool is_non_stator) {
    std::cout << "-- " << path << "\n";
    SearchConfig cfg = parse_config(path);
    Grid grid = build_grid(cfg);

    if (is_non_stator)
        assert(grid.non_stator_positions.count({stator_x, stator_y}) == 1);
    else
        assert(grid.stator_positions.count({stator_x, stator_y}) == 1);

    CadicalSolver solver;
    EncodingStats stats;
    RuleEncoding rule_enc = compute_rule_encoding(cfg.rule);
    TemporalVars tv = encode_all(solver, grid, cfg, stats, rule_enc);
    (void)tv;

    SolverResult result = solver.solve();
    assert(result.status == SolverStatus::SAT);

    // Stator/non-stator must appear alive in the extracted t=0 solution.
    // For stators, encoding must never let them flip dead at any timestep,
    // but extract_full_solution only inspects t=0 — the at-every-t guarantee
    // is already covered by test_stator's forward-sim test. Here we just
    // confirm encoding + extraction agree that the cell is alive in the
    // emitted solution.
    FullSolution sol = extract_full_solution(grid, result, cfg);
    // Bounding box origin is the min of active/unknown/stator/non_stator coords.
    // Recompute it the same way to translate world coords into solution coords.
    int min_x = cfg.pattern.width, min_y = cfg.pattern.height;
    for (int y = 0; y < cfg.pattern.height; y++)
        for (int x = 0; x < cfg.pattern.width; x++) {
            CellState s = cfg.pattern.grid[y][x];
            if (s == CellState::ACTIVE || s == CellState::UNKNOWN ||
                s == CellState::STATOR || s == CellState::NON_STATOR ||
                s == CellState::INIT_OFF) {
                min_x = std::min(min_x, x);
                min_y = std::min(min_y, y);
            }
        }
    assert(sol.cells[stator_y - min_y][stator_x - min_x] == true);
}

int main() {
    std::cout << "=== Stator encoding integration tests ===\n";
    // Both test inputs put the special cell at pattern coord (4, 5).
    run_case("inputs/glider_stator_test.toml",     4, 5, /*non_stator=*/false);
    run_case("inputs/glider_nonstator_test.toml",  4, 5, /*non_stator=*/true);
    std::cout << "=== All stator encoding tests passed ===\n";
    return 0;
}
