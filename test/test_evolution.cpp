#include <iostream>
#include <cassert>
#include <vector>
#include <array>
#include <cstdlib>
#include <string>
#include "sat_evolution.hpp"

// Compute CGOL successor of a grid
std::vector<std::string> cgol_step(const std::vector<std::string>& grid) {
    int h = grid.size();
    int w = grid[0].size();
    std::vector<std::string> next(h, std::string(w, '.'));
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int neighbors = 0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int nx = x + dx, ny = y + dy;
                    if (nx >= 0 && nx < w && ny >= 0 && ny < h && grid[ny][nx] == 'O')
                        neighbors++;
                }
            }
            bool alive = grid[y][x] == 'O';
            if (alive && (neighbors == 2 || neighbors == 3))
                next[y][x] = 'O';
            else if (!alive && neighbors == 3)
                next[y][x] = 'O';
        }
    }
    return next;
}

// Check if all evolution clauses are satisfied given a variable assignment
bool check_evolution_clauses(const PrimeClauseList& clauses,
                              const std::vector<bool>& assignment) {
    for (const auto& clause : clauses) {
        bool sat = false;
        for (int lit : clause) {
            if (lit == 0) continue;
            int var = std::abs(lit);
            bool val = assignment[var - 1];
            if (lit > 0 ? val : !val) {
                sat = true;
                break;
            }
        }
        if (!sat) return false;
    }
    return true;
}

// Test: generate evolution clauses for a 2-gen grid, assign correct successor,
// verify all clauses satisfied. Then flip one cell and verify violation.
void test_evolution(const std::vector<std::string>& gen0, const std::string& name) {
    std::vector<std::string> gen1 = cgol_step(gen0);
    int h = gen0.size();
    int w = gen0[0].size();

    // Assign variable indices: all cells are SAT variables (starting at 2)
    // gen0 cells: var_index = (y * w + x) + 2
    // gen1 cells: var_index = (h * w) + (y * w + x) + 2
    auto var_idx = [&](int x, int y, int t) -> int {
        if (x < 0 || x >= w || y < 0 || y >= h) return 0;  // out-of-bounds = dead
        return t * h * w + y * w + x + 2;
    };

    int total_vars = 2 * h * w;
    std::vector<bool> assignment(total_vars + 1, false);  // 0-indexed SAT vars

    // Set assignment from grids
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            assignment[var_idx(x, y, 0) - 1] = (gen0[y][x] == 'O');
            assignment[var_idx(x, y, 1) - 1] = (gen1[y][x] == 'O');
        }
    }

    // Generate evolution clauses for each cell transition
    PrimeClauseList all_clauses;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            std::array<int, 10> ten;
            int i = 0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    ten[i++] = var_idx(x + dx, y + dy, 0);
                }
            }
            ten[9] = var_idx(x, y, 1);
            auto clauses = generate_evolution_clauses(ten);
            all_clauses.insert(all_clauses.end(), clauses.begin(), clauses.end());
        }
    }

    // Test correct successor
    bool ok = check_evolution_clauses(all_clauses, assignment);
    std::cout << name << " (correct successor): " << (ok ? "PASS" : "FAIL") << "\n";
    assert(ok);

    // Test wrong successor: flip a cell in gen1 that changes
    // Find a cell that differs between gen0 and gen1, or just flip any interior cell
    bool found_violation = false;
    for (int y = 1; y < h - 1 && !found_violation; y++) {
        for (int x = 1; x < w - 1 && !found_violation; x++) {
            int var = var_idx(x, y, 1) - 1;
            assignment[var] = !assignment[var];  // flip

            bool still_ok = check_evolution_clauses(all_clauses, assignment);
            if (!still_ok) {
                std::cout << name << " (flipped cell at " << x << "," << y
                          << "): PASS (violation detected)\n";
                found_violation = true;
            }
            assignment[var] = !assignment[var];  // restore
        }
    }
    assert(found_violation);
}

void test_prime_implicant_count() {
    std::cout << "=== Evolution prime implicant stats ===\n";
    std::cout << "Number of prime implicants: " << evolutionPrimeImplicants.size() << "\n";
    std::cout << "Truth table verification: PASS (checked during static init)\n";
}

int main() {
    test_prime_implicant_count();

    std::cout << "\n=== Testing evolution encoding ===\n";

    // Block (stays the same)
    test_evolution({
        "....",
        ".OO.",
        ".OO.",
        "....",
    }, "Block");

    // Blinker (oscillates)
    test_evolution({
        ".....",
        ".....",
        ".OOO.",
        ".....",
        ".....",
    }, "Blinker");

    // Glider gen 0 → 1
    test_evolution({
        ".....",
        "..O..",
        "...O.",
        ".OOO.",
        ".....",
    }, "Glider");

    std::cout << "\n=== All evolution tests passed ===\n";
    return 0;
}
