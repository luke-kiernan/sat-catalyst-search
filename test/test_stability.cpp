#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include "sat_stability.hpp"
#include "rule.hpp"

// Implicants for the default CGoL rule, used throughout this file.
static const std::vector<std::pair<int,int>>& cgol_implicants() {
    static const auto v = compute_stability_implicants(parse_rule("B3/S23"));
    return v;
}

// Check if all clauses are satisfied given a mapping: SAT_var -> bool
bool all_clauses_satisfied(const StabilityClauseList& clauses,
                           const std::vector<bool>& assignment) {
    for (const auto& clause : clauses) {
        bool sat = false;
        for (int lit : clause) {
            if (lit == 0) continue;
            int var = std::abs(lit);
            bool val = assignment[var - 1];  // var is 1-indexed
            if (lit > 0 ? val : !val) {
                sat = true;
                break;
            }
        }
        if (!sat) return false;
    }
    return true;
}

// Test helper: given a grid (vector of strings, 'O'=alive, '.'=dead),
// check stability for the cell at (cx,cy) using its 3x3 neighborhood.
// All cells get fresh SAT variables (starting from 2).
// Returns whether all stability clauses are satisfied.
bool test_pattern_stability(const std::vector<std::string>& grid, int cx, int cy,
                            const std::string& name) {
    int h = grid.size();
    int w = grid[0].size();

    // Assign variable indices: 0=dead outside, otherwise fresh var starting at 2
    auto cell_index = [&](int x, int y) -> int {
        if (x < 0 || x >= w || y < 0 || y >= h) return 0;  // dead
        return (y * w + x) + 2;  // SAT var starts at 2
    };

    // Build 9-cell neighborhood for (cx, cy)
    std::array<int, 9> nine;
    int i = 0;
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            nine[i++] = cell_index(cx + dx, cy + dy);
        }
    }

    StabilityClauseList clauses = generate_stability_clauses(nine, cgol_implicants());

    // Build assignment from grid
    int max_var = h * w + 1;  // max SAT var index
    std::vector<bool> assignment(max_var, false);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int var = cell_index(x, y) - 1;  // 0-indexed SAT var
            assignment[var] = (grid[y][x] == 'O');
        }
    }

    return all_clauses_satisfied(clauses, assignment);
}

// Test full pattern: every cell in the interior must be stable
bool test_full_pattern_stability(const std::vector<std::string>& grid,
                                  const std::string& name) {
    int h = grid.size();
    int w = grid[0].size();

    for (int cy = 0; cy < h; cy++) {
        for (int cx = 0; cx < w; cx++) {
            if (!test_pattern_stability(grid, cx, cy, name)) {
                std::cout << "  UNSTABLE at (" << cx << "," << cy << ")\n";
                return false;
            }
        }
    }
    return true;
}

void test_still_lifes() {
    std::cout << "=== Testing still lifes (should all be stable) ===\n";

    // Block (2x2)
    {
        std::vector<std::string> grid = {
            "....",
            ".OO.",
            ".OO.",
            "....",
        };
        bool ok = test_full_pattern_stability(grid, "block");
        std::cout << "Block: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }

    // Beehive
    {
        std::vector<std::string> grid = {
            "......",
            "..OO..",
            ".O..O.",
            "..OO..",
            "......",
        };
        bool ok = test_full_pattern_stability(grid, "beehive");
        std::cout << "Beehive: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }

    // Loaf
    {
        std::vector<std::string> grid = {
            "......",
            "..OO..",
            ".O..O.",
            "..O.O.",
            "...O..",
            "......",
        };
        bool ok = test_full_pattern_stability(grid, "loaf");
        std::cout << "Loaf: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }

    // Boat
    {
        std::vector<std::string> grid = {
            ".....",
            ".OO..",
            ".O.O.",
            "..O..",
            ".....",
        };
        bool ok = test_full_pattern_stability(grid, "boat");
        std::cout << "Boat: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }

    // Tub
    {
        std::vector<std::string> grid = {
            ".....",
            "..O..",
            ".O.O.",
            "..O..",
            ".....",
        };
        bool ok = test_full_pattern_stability(grid, "tub");
        std::cout << "Tub: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }

    // All-S23-transitions still life (exercises every survival case)
    {
        std::vector<std::string> grid = {
            ".............",
            "...OO........",
            "....O........",
            "...O.....OO..",
            "..O...O..OO..",
            ".O.OOOO......",
            "..O....OOOO..",
            "...OOO..O..O.",
            ".....OO...OO.",
            ".............",
        };
        bool ok = test_full_pattern_stability(grid, "all-S23");
        std::cout << "All-S23-transitions: " << (ok ? "PASS" : "FAIL") << "\n";
        assert(ok);
    }
}

void test_non_still_lifes() {
    std::cout << "\n=== Testing non-still-lifes (should be unstable) ===\n";

    // Blinker (horizontal phase)
    {
        std::vector<std::string> grid = {
            ".....",
            ".....",
            ".OOO.",
            ".....",
            ".....",
        };
        bool ok = test_full_pattern_stability(grid, "blinker");
        std::cout << "Blinker: " << (!ok ? "PASS" : "FAIL") << " (should be unstable)\n";
        assert(!ok);
    }

    // Glider
    {
        std::vector<std::string> grid = {
            ".....",
            "..O..",
            "...O.",
            ".OOO.",
            ".....",
        };
        bool ok = test_full_pattern_stability(grid, "glider");
        std::cout << "Glider: " << (!ok ? "PASS" : "FAIL") << " (should be unstable)\n";
        assert(!ok);
    }

    // R-pentomino
    {
        std::vector<std::string> grid = {
            ".....",
            "..OO.",
            ".OO..",
            "..O..",
            ".....",
        };
        bool ok = test_full_pattern_stability(grid, "R-pentomino");
        std::cout << "R-pentomino: " << (!ok ? "PASS" : "FAIL") << " (should be unstable)\n";
        assert(!ok);
    }

    // Isolated cell (1 alive cell — dies from underpopulation)
    {
        std::vector<std::string> grid = {
            "...",
            ".O.",
            "...",
        };
        bool ok = test_full_pattern_stability(grid, "isolated cell");
        std::cout << "Isolated cell: " << (!ok ? "PASS" : "FAIL") << " (should be unstable)\n";
        assert(!ok);
    }

    // Block with one cell removed (3 cells in L — unstable)
    {
        std::vector<std::string> grid = {
            "....",
            ".OO.",
            ".O..",
            "....",
        };
        bool ok = test_full_pattern_stability(grid, "block-minus-one");
        std::cout << "Block minus one cell: " << (!ok ? "PASS" : "FAIL") << " (should be unstable)\n";
        assert(!ok);
    }
}

void test_prime_implicant_count() {
    std::cout << "\n=== Stability prime implicant stats ===\n";
    std::cout << "Number of prime implicants: " << cgol_implicants().size() << "\n";
    // The assertion in the static initializer already verified correctness for all 512 configs
    std::cout << "Truth table verification: PASS (checked during static init)\n";
}

int main() {
    test_prime_implicant_count();
    test_still_lifes();
    test_non_still_lifes();
    std::cout << "\n=== All stability tests passed ===\n";
    return 0;
}
