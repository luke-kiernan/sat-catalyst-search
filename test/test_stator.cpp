#include <iostream>
#include <cassert>
#include <string>
#include "grid.hpp"
#include "config.hpp"

// Build an in-memory SearchConfig from an inline RLE — avoids per-test TOML files.
static SearchConfig make_config(const std::string& rle) {
    SearchConfig cfg;
    cfg.pattern = parse_rle(rle);
    cfg.center_x = cfg.pattern.width / 2;
    cfg.center_y = cfg.pattern.height / 2;
    cfg.first_active_range = {1, 5};
    cfg.active_window_range = {1, 5};
    cfg.min_stable_interval = 2;
    return cfg;
}

// 1. RLE parser routes C/D/E through the switch into the Barrister-matching
//    states: C (3) = NON_STATOR, D (4) = INIT_OFF, E (5) = STATOR.
static void test_rle_parses_catalyst_states() {
    std::cout << "-- RLE parses C/D/E into NON_STATOR/INIT_OFF/STATOR\n";
    auto p = parse_rle("x = 5, y = 1, rule = LifeHistory\nCBEDB!");
    assert(p.grid[0][0] == CellState::NON_STATOR);
    assert(p.grid[0][1] == CellState::UNKNOWN);
    assert(p.grid[0][2] == CellState::STATOR);
    assert(p.grid[0][3] == CellState::INIT_OFF);
    assert(p.grid[0][4] == CellState::UNKNOWN);
}

// 2. build_grid routes STATOR/NON_STATOR to the known-alive branch:
//    populates the right position sets, and assigns catalyst_var = 1
//    instead of a fresh SAT variable.
static void test_grid_dispatches_stator_and_non_stator() {
    std::cout << "-- build_grid routes E/C to stator/non_stator sets with var=1\n";
    // Pattern: glider rows 0-2, then 3 blank rows, then BEB / BCB at rows 6,7.
    // E (state 5) = stator at (5,6); C (state 3) = non_stator at (5,7).
    SearchConfig cfg = make_config(
        "x = 7, y = 8, rule = LifeHistory\n"
        ".A$2.A$3A$3$4.BEB$4.BCB!");
    Grid g = build_grid(cfg);

    assert(g.stator_positions.count({5, 6}) == 1);
    assert(g.non_stator_positions.count({5, 7}) == 1);
    assert(g.catalyst_positions.count({5, 6}) == 0);  // NOT unknown
    assert(g.catalyst_positions.count({5, 7}) == 0);

    // Known-alive sentinel, not a fresh SAT variable.
    assert(g.catalyst_var_at(5, 6) == 1);
    assert(g.catalyst_var_at(5, 7) == 1);
    // Adjacent unknown cells should get real SAT variables.
    assert(g.catalyst_var_at(4, 6) >= 2);
}

// 3. Headline behavior: stator cells are pinned alive at every t in the
//    forward simulator, even when the neighborhood alone wouldn't keep them alive.
static void test_stator_pinned_alive_across_timesteps() {
    std::cout << "-- Stator stays alive at every t even when neighbors wouldn't sustain it\n";
    // Single stator with only B cells around it. Without the pin, the forward
    // simulator would mark it unknown/dead at t>=1 because evolution depends on
    // the unknown neighbors. With the pin, it stays =1 at every timestep.
    // Glider rows 0-2, blank row, then 3B / BEB / 3B at rows 3,4,5. Stator (E) at (5,4).
    SearchConfig cfg = make_config(
        "x = 8, y = 6, rule = LifeHistory\n"
        ".A$2.A$3A$4.3B$4.BEB$4.3B!");
    Grid g = build_grid(cfg);

    for (int t = 0; t < g.total_gens; t++) {
        assert(g.cell_at(5, 4, t) == 1);
    }
}

// TODO(Phase 5): The "isolated stator → throw at build_grid" check was driven
// by the old StableMaskGrid forward sim, which Phase 4 removed. The same
// situation now manifests as UNSAT at solve time (the evolution constraint
// `1 = rule.evolves_to(stator=1, 8 dead)` is false under B3/S23). When we
// rebuild the precomputation in a way that re-introduces forward analysis,
// or add an explicit pre-encoding sanity check, restore this assertion.
static void test_isolated_stator_throws() {
    std::cout << "-- (skipped) Isolated stator throws — see TODO\n";
}

// 5. ZOI computation unions all three catalyst sets — perturbation_region
//    and catalyst_neighborhood must cover stator/non_stator cells and their
//    8-neighborhoods. Regression guard for "forgot to add a new cell type to
//    all_catalyst" — without this, recovery/perturbation constraints would
//    silently skip the stator's surroundings.
static void test_zoi_includes_stator_and_non_stator() {
    std::cout << "-- Perturbation region/neighborhood cover stator+non_stator ZOI\n";
    // Same layout as test 2: stator (E) at (5,6), non_stator (C) at (5,7).
    SearchConfig cfg = make_config(
        "x = 7, y = 8, rule = LifeHistory\n"
        ".A$2.A$3A$3$4.BEB$4.BCB!");
    Grid g = build_grid(cfg);

    // Diagonals of the stator/non_stator cells — only included if the ZOI
    // loop iterated over stator_positions and non_stator_positions.
    assert(g.perturbation_region.count({4, 5}) == 1);  // stator NW diagonal
    assert(g.perturbation_region.count({6, 5}) == 1);  // stator NE diagonal
    assert(g.perturbation_region.count({4, 8}) == 1);  // non_stator SW diagonal
    assert(g.perturbation_region.count({6, 8}) == 1);  // non_stator SE diagonal

    // catalyst_neighborhood contains the cells themselves
    assert(g.catalyst_neighborhood.count({5, 6}) == 1);
    assert(g.catalyst_neighborhood.count({5, 7}) == 1);
}

// 6. Headline behavior for INIT_OFF (state 4 / D): the cell is a known
//    catalyst cell (stable value alive, var=1, in the ZOI) but its gen-0
//    state is DEAD — the distinguishing property versus NON_STATOR, which is
//    alive at gen 0. This is the behavior that lets a catalyst "fill in".
static void test_init_off_dead_at_gen0_but_known() {
    std::cout << "-- INIT_OFF is a known catalyst cell but dead at gen 0\n";
    // Glider rows 0-2, blank row, then 3B / BDB / 3B at rows 3,4,5. INIT_OFF (D) at (5,4).
    SearchConfig cfg = make_config(
        "x = 8, y = 6, rule = LifeHistory\n"
        ".A$2.A$3A$4.3B$4.BDB$4.3B!");
    Grid g = build_grid(cfg);

    // Routed to init_off, not the other catalyst sets.
    assert(g.init_off_positions.count({5, 4}) == 1);
    assert(g.stator_positions.count({5, 4}) == 0);
    assert(g.non_stator_positions.count({5, 4}) == 0);
    assert(g.catalyst_positions.count({5, 4}) == 0);  // not an unknown cell

    // Stable value is alive (known-alive sentinel), so it conditions ZOI.
    assert(g.catalyst_var_at(5, 4) == 1);
    assert(g.catalyst_neighborhood.count({5, 4}) == 1);

    // The distinguishing behavior: dead at generation 0.
    assert(g.cell_at(5, 4, 0) == 0);
}

int main() {
    std::cout << "=== Stator / non-stator tests ===\n";
    test_rle_parses_catalyst_states();
    test_grid_dispatches_stator_and_non_stator();
    test_stator_pinned_alive_across_timesteps();
    test_isolated_stator_throws();
    test_zoi_includes_stator_and_non_stator();
    test_init_off_dead_at_gen0_but_known();
    std::cout << "=== All stator tests passed ===\n";
    return 0;
}
