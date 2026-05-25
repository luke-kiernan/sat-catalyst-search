#include <iostream>
#include <cassert>
#include "config.hpp"

int main() {
    std::cout << "=== Testing RLE parser + TOML config ===\n\n";

    SearchConfig config = parse_config("inputs/glider.toml");

    std::cout << "Pattern dimensions: " << config.pattern.width << " x "
              << config.pattern.height << "\n";
    std::cout << "Pattern center: (" << config.center_x << ", " << config.center_y << ")\n";
    std::cout << "First active range: [" << config.first_active_range[0] << ", "
              << config.first_active_range[1] << "]\n";
    std::cout << "Active window range: [" << config.active_window_range[0] << ", "
              << config.active_window_range[1] << "]\n";
    std::cout << "Min stable interval: " << config.min_stable_interval << "\n";
    std::cout << "Ever active bounds: [" << config.ever_active_bounds[0] << ", "
              << config.ever_active_bounds[1] << "]\n";
    std::cout << "Total generations: " << config.total_generations() << "\n\n";

    // Count cells by type
    int dead = 0, active = 0, unknown = 0;
    for (int y = 0; y < config.pattern.height; y++) {
        for (int x = 0; x < config.pattern.width; x++) {
            switch (config.pattern.grid[y][x]) {
                case CellState::DEAD: dead++; break;
                case CellState::ACTIVE: active++; break;
                case CellState::UNKNOWN: unknown++; break;
                default: break;
            }
        }
    }
    std::cout << "Cell counts: dead=" << dead << " active=" << active
              << " unknown=" << unknown << "\n";

    // Verify against manual inspection of glider.toml:
    // 29x24 grid, glider at top-left, big rectangle of B (unknown) cells
    assert(config.pattern.width == 29);
    assert(config.pattern.height == 24);
    assert(config.center_x == 15);
    assert(config.center_y == 10);
    assert(config.first_active_range[0] == 3);
    assert(config.first_active_range[1] == 7);
    assert(config.active_window_range[0] == 3);
    assert(config.active_window_range[1] == 40);
    assert(config.min_stable_interval == 5);
    assert(active == 5);  // glider has 5 cells
    assert(unknown > 0);   // should have many unknown cells

    // Print pattern
    std::cout << "\nPattern grid:\n";
    for (int y = 0; y < config.pattern.height; y++) {
        for (int x = 0; x < config.pattern.width; x++) {
            switch (config.pattern.grid[y][x]) {
                case CellState::DEAD: std::cout << '.'; break;
                case CellState::ACTIVE: std::cout << 'A'; break;
                case CellState::UNKNOWN: std::cout << 'B'; break;
                case CellState::NON_STATOR: std::cout << 'C'; break;
                case CellState::INIT_OFF: std::cout << 'D'; break;
                case CellState::STATOR: std::cout << 'E'; break;
            }
        }
        std::cout << "\n";
    }

    std::cout << "\n=== Config test passed ===\n";
    return 0;
}
