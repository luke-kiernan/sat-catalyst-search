#include <iostream>
#include <vector>
#include <string>

std::vector<std::string> step(const std::vector<std::string>& grid) {
    int h = grid.size(), w = grid[0].size();
    std::vector<std::string> next(h, std::string(w, '.'));
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            int n = 0;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;
                    int nx = x+dx, ny = y+dy;
                    if (nx>=0&&nx<w&&ny>=0&&ny<h&&grid[ny][nx]=='O') n++;
                }
            bool alive = grid[y][x] == 'O';
            if ((alive && (n==2||n==3)) || (!alive && n==3))
                next[y][x] = 'O';
        }
    return next;
}

void show(const std::string& label, const std::vector<std::string>& g) {
    std::cout << label << ":\n";
    for (auto& r : g) std::cout << "  " << r << "\n";
    std::cout << "\n";
}

int main() {
    // Pattern from barrister glider.toml: .O./..O/OOO
    std::vector<std::string> g1 = {
        "......",
        "..O...",
        "...O..",
        ".OOO..",
        "......",
        "......",
    };
    show("barrister pattern gen 0", g1);
    for (int i = 1; i <= 4; i++) {
        g1 = step(g1);
        show("gen " + std::to_string(i), g1);
    }

    // Canonical glider: .O/O./OOO
    std::cout << "---\n\n";
    std::vector<std::string> g2 = {
        "......",
        "..O...",
        ".O....",
        ".OOO..",
        "......",
        "......",
    };
    show("canonical bo$o$3o gen 0", g2);
    for (int i = 1; i <= 4; i++) {
        g2 = step(g2);
        show("gen " + std::to_string(i), g2);
    }
}
