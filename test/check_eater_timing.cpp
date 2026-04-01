#include <iostream>
#include <vector>
#include <string>

using Grid = std::vector<std::string>;

Grid step(const Grid& g) {
    int h = g.size(), w = g[0].size();
    Grid next(h, std::string(w, '.'));
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            int n = 0;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int nx = x+dx, ny = y+dy;
                    if (nx>=0&&nx<w&&ny>=0&&ny<h&&g[ny][nx]=='O') n++;
                }
            bool alive = g[y][x] == 'O';
            if ((alive && (n==2||n==3)) || (!alive && n==3))
                next[y][x] = 'O';
        }
    return next;
}

void show(const Grid& g, int gen) {
    std::cout << "Gen " << gen << ":\n";
    for (auto& r : g) {
        for (char c : r) std::cout << (c == 'O' ? 'O' : '.');
        std::cout << "\n";
    }
}

int main() {
    int W = 15, H = 15;

    // Glider: .O / ..O / OOO at (1,0),(2,1),(0,2),(1,2),(2,2)
    auto make_grid = [&]() -> Grid {
        Grid g(H, std::string(W, '.'));
        g[0][1] = 'O';
        g[1][2] = 'O';
        g[2][0] = 'O'; g[2][1] = 'O'; g[2][2] = 'O';
        return g;
    };

    // Eater-1: OO / O.O / ..O / ..OO
    auto place_eater = [](Grid& g, int ex, int ey) {
        g[ey][ex] = 'O'; g[ey][ex+1] = 'O';
        g[ey+1][ex] = 'O'; g[ey+1][ex+2] = 'O';
        g[ey+2][ex+2] = 'O';
        g[ey+3][ex+2] = 'O'; g[ey+3][ex+3] = 'O';
    };

    // Perturbation region: ZOI of catalyst cells (4,4)-(8,8)
    // = cells (3,3)-(9,9)
    auto is_catalyst = [](int x, int y) {
        return x >= 4 && x <= 8 && y >= 4 && y <= 8;
    };
    auto in_perturbation = [](int x, int y) {
        return x >= 3 && x <= 9 && y >= 3 && y <= 9;
    };

    // Stable eater assignments (which catalyst cells are alive)
    auto eater44_stable = [](int x, int y) -> bool {
        // Eater at (4,4): cells (4,4),(5,4),(4,5),(6,5),(6,6),(6,7),(7,7)
        return (x==4&&y==4)||(x==5&&y==4)||(x==4&&y==5)||(x==6&&y==5)||
               (x==6&&y==6)||(x==6&&y==7)||(x==7&&y==7);
    };
    auto eater55_stable = [](int x, int y) -> bool {
        // Eater at (5,5): cells (5,5),(6,5),(5,6),(7,6),(7,7),(7,8),(8,8)
        return (x==5&&y==5)||(x==6&&y==5)||(x==5&&y==6)||(x==7&&y==6)||
               (x==7&&y==7)||(x==7&&y==8)||(x==8&&y==8);
    };

    // Free evolution (no catalyst)
    Grid free_grid = make_grid();

    // Eater at (4,4)
    Grid eater44 = make_grid();
    place_eater(eater44, 4, 4);

    // Eater at (5,5)
    Grid eater55 = make_grid();
    place_eater(eater55, 5, 5);

    std::cout << "=== Comparing free evolution vs eater@(4,4) vs eater@(5,5) ===\n";
    std::cout << "Perturbation region: (3,3)-(9,9), catalyst region: (4,4)-(8,8)\n";
    std::cout << "Catalyst cells: compare against stable eater assignment\n";
    std::cout << "Non-catalyst cells: compare against free evolution\n\n";

    for (int gen = 0; gen <= 12; gen++) {
        // Check perturbation using SAT-matching logic
        bool diff44 = false, diff55 = false;
        std::string where44, where55;
        for (int y = 3; y <= 9; y++) {
            for (int x = 3; x <= 9; x++) {
                if (is_catalyst(x, y)) {
                    // Catalyst cell: compare against stable eater assignment
                    bool actual44 = (eater44[y][x] == 'O');
                    bool stable44 = eater44_stable(x, y);
                    if (actual44 != stable44) {
                        diff44 = true;
                        where44 += " cat(" + std::to_string(x) + "," + std::to_string(y) + ")";
                    }
                    bool actual55 = (eater55[y][x] == 'O');
                    bool stable55 = eater55_stable(x, y);
                    if (actual55 != stable55) {
                        diff55 = true;
                        where55 += " cat(" + std::to_string(x) + "," + std::to_string(y) + ")";
                    }
                } else {
                    // Non-catalyst cell: compare against free evolution
                    if (free_grid[y][x] != eater44[y][x]) {
                        diff44 = true;
                        where44 += " free(" + std::to_string(x) + "," + std::to_string(y) + ")";
                    }
                    if (free_grid[y][x] != eater55[y][x]) {
                        diff55 = true;
                        where55 += " free(" + std::to_string(x) + "," + std::to_string(y) + ")";
                    }
                }
            }
        }

        std::cout << "Gen " << gen << ":";
        std::cout << "  @(4,4)=" << (diff44 ? "YES" : "no ");
        std::cout << "  @(5,5)=" << (diff55 ? "YES" : "no ");
        if (diff44) std::cout << " [" << where44 << " ]";
        if (diff55) std::cout << " [" << where55 << " ]";
        std::cout << "\n";

        if (diff44 || diff55) {
            std::cout << "  Free:          Eater@(4,4):     Eater@(5,5):\n";
            for (int y = 0; y < H; y++) {
                std::cout << "  ";
                for (int x = 0; x < W; x++)
                    std::cout << (free_grid[y][x] == 'O' ? 'O' : '.');
                std::cout << "  ";
                for (int x = 0; x < W; x++)
                    std::cout << (eater44[y][x] == 'O' ? 'O' : '.');
                std::cout << "  ";
                for (int x = 0; x < W; x++)
                    std::cout << (eater55[y][x] == 'O' ? 'O' : '.');
                std::cout << "\n";
            }
            std::cout << "\n";
        }

        free_grid = step(free_grid);
        eater44 = step(eater44);
        eater55 = step(eater55);
    }
}
