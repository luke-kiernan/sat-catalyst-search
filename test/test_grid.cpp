#include <iostream>
#include <cassert>
#include "grid.hpp"

int main() {
    std::cout << "=== Testing grid construction ===\n\n";

    SearchConfig config = parse_config("inputs/glider.toml");
    Grid grid = build_grid(config);

    std::cout << "Total generations: " << grid.total_gens << "\n";
    std::cout << "Grid dimensions: " << grid.width << " x " << grid.height << "\n";
    std::cout << "Grid origin: (" << grid.ox << ", " << grid.oy << ")\n";
    std::cout << "Total SAT variables: " << grid.num_vars() << "\n";
    std::cout << "Catalyst positions: " << grid.catalyst_positions.size() << "\n";
    std::cout << "Perturbation region: " << grid.perturbation_region.size() << "\n";
    std::cout << "Catalyst neighborhood: " << grid.catalyst_neighborhood.size() << "\n";

    // Verify t=0 has correct known cells
    // Glider at pattern coordinates: (7,1)=A, (8,2)=A, (6,3)=A, (7,3)=A, (8,3)=A
    assert(grid.cell_at(7, 1, 0) == 1);  // active
    assert(grid.cell_at(8, 2, 0) == 1);  // active
    assert(grid.cell_at(6, 3, 0) == 1);  // active
    assert(grid.cell_at(7, 3, 0) == 1);  // active
    assert(grid.cell_at(8, 3, 0) == 1);  // active
    std::cout << "\nGlider cells at t=0: PASS (all alive)\n";

    // Verify dead cells
    assert(grid.cell_at(0, 0, 0) == 0);  // dead
    std::cout << "Dead cell at (0,0,0): PASS\n";

    // Verify catalyst cells are SAT variables (>=2)
    assert(grid.cell_at(11, 0, 0) >= 2);  // first unknown column
    std::cout << "Catalyst cell at (11,0,0): var " << grid.cell_at(11, 0, 0) << " PASS\n";

    // Verify light cone expansion: count non-zero cells per timestep
    std::cout << "\nLight cone sizes (non-zero cells per timestep):\n";
    for (int t = 0; t < std::min(5, grid.total_gens); t++) {
        int count = 0;
        for (int y = 0; y < grid.height; y++) {
            for (int x = 0; x < grid.width; x++) {
                if (grid.cells[t][y][x] != 0) count++;
            }
        }
        std::cout << "  t=" << t << ": " << count << " cells\n";
    }

    // Light cone should expand each generation
    auto count_cells = [&](int t) {
        int c = 0;
        for (int y = 0; y < grid.height; y++)
            for (int x = 0; x < grid.width; x++)
                if (grid.cells[t][y][x] != 0) c++;
        return c;
    };

    for (int t = 1; t < std::min(5, grid.total_gens); t++) {
        assert(count_cells(t) >= count_cells(t - 1));
    }
    std::cout << "Light cone expansion: PASS\n";

    std::cout << "\n=== Grid test passed ===\n";
    return 0;
}
