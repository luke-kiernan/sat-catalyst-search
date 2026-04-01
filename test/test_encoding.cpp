#include <iostream>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <set>
#include "encoding.hpp"
#include "output.hpp"

// Compute CGOL successor
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

// Simple test: build the grid and encoding for a known pattern,
// then manually check if the clauses are satisfiable with a known catalyst.
int main() {
    std::cout << "=== Encoding diagnostic test ===\n\n";

    // Use the tiny glider input
    SearchConfig config = parse_config("inputs/glider_tiny.toml");

    std::cout << "Config loaded:\n";
    std::cout << "  Pattern: " << config.pattern.width << "x" << config.pattern.height << "\n";
    std::cout << "  first-active-range: [" << config.first_active_range[0] << ", "
              << config.first_active_range[1] << "]\n";
    std::cout << "  active-window-range: [" << config.active_window_range[0] << ", "
              << config.active_window_range[1] << "]\n";
    std::cout << "  min-stable-interval: " << config.min_stable_interval << "\n";
    std::cout << "  total_generations: " << config.total_generations() << "\n\n";

    // Print the parsed pattern
    std::cout << "Parsed pattern:\n";
    for (int y = 0; y < config.pattern.height; y++) {
        for (int x = 0; x < config.pattern.width; x++) {
            switch (config.pattern.grid[y][x]) {
                case CellState::DEAD: std::cout << '.'; break;
                case CellState::ACTIVE: std::cout << 'A'; break;
                case CellState::UNKNOWN: std::cout << 'B'; break;
                default: std::cout << '?'; break;
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Now simulate the glider + a block catalyst at (5,2) to see what happens
    // Block at columns 5-6, rows 2-3
    int sim_w = 15, sim_h = 12;
    std::vector<std::string> sim(sim_h, std::string(sim_w, '.'));

    // Place glider (from pattern: A at (1,0), A at (2,1), AAA at (0,2)(1,2)(2,2))
    sim[0][1] = 'O'; // (1,0)
    sim[1][2] = 'O'; // (2,1)
    sim[2][0] = 'O'; // (0,2)
    sim[2][1] = 'O'; // (1,2)
    sim[2][2] = 'O'; // (2,2)

    // Place a block catalyst at (5,1)-(6,2)
    sim[1][5] = 'O';
    sim[1][6] = 'O';
    sim[2][5] = 'O';
    sim[2][6] = 'O';

    std::cout << "Simulation (glider + block):\n";
    int total_gens = config.total_generations() + 1;
    for (int t = 0; t < std::min(total_gens, 20); t++) {
        std::cout << "Gen " << t << ":\n";
        for (int y = 0; y < sim_h; y++) {
            for (int x = 0; x < sim_w; x++) {
                std::cout << (sim[y][x] == 'O' ? 'O' : '.');
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        // Check if block at (5,1)-(6,2) is intact
        bool block_ok = sim[1][5] == 'O' && sim[1][6] == 'O' &&
                        sim[2][5] == 'O' && sim[2][6] == 'O';
        // Check neighbors of block are dead
        bool neighbors_clear = true;
        for (int by = 0; by <= 3; by++) {
            for (int bx = 4; bx <= 7; bx++) {
                if ((bx >= 5 && bx <= 6 && by >= 1 && by <= 2)) continue; // block itself
                if (sim[by][bx] == 'O') neighbors_clear = false;
            }
        }

        std::cout << "  Block intact: " << (block_ok ? "yes" : "NO")
                  << ", neighbors clear: " << (neighbors_clear ? "yes" : "NO") << "\n\n";

        sim = cgol_step(sim);
    }

    std::cout << "=== Diagnostic complete ===\n";
    return 0;
}
