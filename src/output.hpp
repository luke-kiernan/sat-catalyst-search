#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "grid.hpp"
#include "config.hpp"
#include "solver_result.hpp"

// ── Summary RLE output ──
// Tile all solutions into a single RLE grid, arranged in columns of SUMMARY_ROWS
// with spacing aligned to multiples of 10.

constexpr int SUMMARY_COLS = 8;
constexpr int SUMMARY_CELL_BUFFER = 10;

// Extract catalyst pattern from solver solution
struct CatalystSolution {
    std::vector<std::vector<bool>> cells;  // [y][x] relative to catalyst bbox
    int min_x, min_y, max_x, max_y;       // world coordinates of bounding box
    int population = 0;
};

inline CatalystSolution extract_catalyst(const Grid& grid, const SolverResult& result) {
    CatalystSolution sol;
    sol.min_x = 1 << 30;
    sol.min_y = 1 << 30;
    sol.max_x = -(1 << 30);
    sol.max_y = -(1 << 30);
    sol.population = 0;

    // First pass: find bounding box of alive catalyst cells
    for (auto [wx, wy] : grid.catalyst_positions) {
        int var = grid.catalyst_var_at(wx, wy);
        if (var < 2) continue;
        int sat_var = var - 1;
        bool alive = result.solution.count(sat_var) > 0;
        if (alive) {
            sol.min_x = std::min(sol.min_x, wx);
            sol.min_y = std::min(sol.min_y, wy);
            sol.max_x = std::max(sol.max_x, wx);
            sol.max_y = std::max(sol.max_y, wy);
            sol.population++;
        }
    }

    if (sol.population == 0) return sol;

    // Second pass: fill grid
    int w = sol.max_x - sol.min_x + 1;
    int h = sol.max_y - sol.min_y + 1;
    sol.cells.resize(h, std::vector<bool>(w, false));

    for (auto [wx, wy] : grid.catalyst_positions) {
        int var = grid.catalyst_var_at(wx, wy);
        if (var < 2) continue;
        int sat_var = var - 1;
        bool alive = result.solution.count(sat_var) > 0;
        if (alive) {
            sol.cells[wy - sol.min_y][wx - sol.min_x] = true;
        }
    }

    return sol;
}

// Convert catalyst solution to RLE string
inline std::string catalyst_to_rle(const CatalystSolution& sol) {
    if (sol.population == 0) return "";

    int h = sol.cells.size();
    int w = sol.cells[0].size();

    std::string rle = "x = " + std::to_string(w) + ", y = " + std::to_string(h)
                    + ", rule = B3/S23\n";

    for (int y = 0; y < h; y++) {
        if (y > 0) rle += "$";

        int x = 0;
        while (x < w) {
            bool state = sol.cells[y][x];
            int run = 1;
            while (x + run < w && sol.cells[y][x + run] == state) run++;

            char c = state ? 'o' : 'b';
            // Skip trailing dead cells
            if (!state && x + run >= w) break;

            if (run > 1) rle += std::to_string(run);
            rle += c;
            x += run;
        }
    }
    rle += "!";

    return rle;
}

// Print catalyst as visual grid
inline void print_catalyst(const CatalystSolution& sol) {
    if (sol.population == 0) {
        std::cout << "(empty)\n";
        return;
    }
    for (const auto& row : sol.cells) {
        for (bool c : row) {
            std::cout << (c ? 'O' : '.');
        }
        std::cout << "\n";
    }
}

// ── Full solution extraction (active pattern + catalyst) ──

struct FullSolution {
    std::vector<std::vector<bool>> cells; // [y][x]
    std::vector<std::vector<bool>> relevant; // [y][x] — true if cell is relevant
    int width = 0, height = 0;
};

// Extract the full t=0 state: active cells + resolved catalyst.
// Bounding box is the tight extent of active + unknown cells in the input pattern.
// If relevant_cells is provided, marks which cells are relevant for blocking.
inline FullSolution extract_full_solution(const Grid& grid, const SolverResult& result,
                                           const SearchConfig& config,
                                           const std::set<std::pair<int,int>>* relevant_cells = nullptr) {
    const auto& pattern = config.pattern;
    int min_x = pattern.width, min_y = pattern.height;
    int max_x = -1, max_y = -1;

    for (int y = 0; y < pattern.height; y++) {
        for (int x = 0; x < pattern.width; x++) {
            CellState s = pattern.grid[y][x];
            if (s == CellState::ACTIVE || s == CellState::UNKNOWN ||
                s == CellState::FORCED_ON || s == CellState::FORCED_ON_INIT_OFF) {
                min_x = std::min(min_x, x);
                min_y = std::min(min_y, y);
                max_x = std::max(max_x, x);
                max_y = std::max(max_y, y);
            }
        }
    }

    FullSolution sol;
    if (max_x < 0) return sol;

    sol.width = max_x - min_x + 1;
    sol.height = max_y - min_y + 1;
    sol.cells.resize(sol.height, std::vector<bool>(sol.width, false));
    sol.relevant.resize(sol.height, std::vector<bool>(sol.width, false));

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            CellState s = pattern.grid[y][x];
            if (s == CellState::ACTIVE || s == CellState::FORCED_ON) {
                sol.cells[y - min_y][x - min_x] = true;
            } else if (s == CellState::UNKNOWN) {
                int var = grid.catalyst_var_at(x, y);
                if (var >= 2)
                    sol.cells[y - min_y][x - min_x] = result.solution.count(var - 1) > 0;
            }
            if (relevant_cells && relevant_cells->count({x, y}))
                sol.relevant[y - min_y][x - min_x] = true;
        }
    }

    return sol;
}

// Extract full solution using a CatalystSolution (e.g. minimized) for catalyst cells.
inline FullSolution extract_full_solution_from_catalyst(const CatalystSolution& cat,
                                                         const SearchConfig& config) {
    const auto& pattern = config.pattern;
    int min_x = pattern.width, min_y = pattern.height;
    int max_x = -1, max_y = -1;

    for (int y = 0; y < pattern.height; y++) {
        for (int x = 0; x < pattern.width; x++) {
            CellState s = pattern.grid[y][x];
            if (s == CellState::ACTIVE || s == CellState::UNKNOWN ||
                s == CellState::FORCED_ON || s == CellState::FORCED_ON_INIT_OFF) {
                min_x = std::min(min_x, x);
                min_y = std::min(min_y, y);
                max_x = std::max(max_x, x);
                max_y = std::max(max_y, y);
            }
        }
    }

    FullSolution sol;
    if (max_x < 0) return sol;

    sol.width = max_x - min_x + 1;
    sol.height = max_y - min_y + 1;
    sol.cells.resize(sol.height, std::vector<bool>(sol.width, false));
    sol.relevant.resize(sol.height, std::vector<bool>(sol.width, false));

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            CellState s = pattern.grid[y][x];
            if (s == CellState::ACTIVE || s == CellState::FORCED_ON) {
                sol.cells[y - min_y][x - min_x] = true;
            } else if (s == CellState::UNKNOWN) {
                // Look up in the CatalystSolution
                if (cat.population > 0 &&
                    x >= cat.min_x && x <= cat.max_x &&
                    y >= cat.min_y && y <= cat.max_y &&
                    cat.cells[y - cat.min_y][x - cat.min_x]) {
                    sol.cells[y - min_y][x - min_x] = true;
                }
            }
        }
    }

    return sol;
}

// Round up to the next multiple of `m`.
inline int round_up_to(int val, int m) {
    return ((val + m - 1) / m) * m;
}

// Write all solutions tiled into a single RLE.
// Layout: SUMMARY_COLS solutions per row, rows grow downward.
// Spacing aligned to multiples of 10.
inline void write_summary_rle(const std::vector<FullSolution>& solutions,
                               std::ostream& out) {
    if (solutions.empty()) return;

    int cell_w = solutions[0].width;
    int cell_h = solutions[0].height;
    int x_spacing = round_up_to(cell_w + SUMMARY_CELL_BUFFER, 10);
    int y_spacing = round_up_to(cell_h + SUMMARY_CELL_BUFFER, 10);

    int n = (int)solutions.size();
    int num_cols = std::min(n, SUMMARY_COLS);
    int num_rows = (n + SUMMARY_COLS - 1) / SUMMARY_COLS;

    int total_w = num_cols * x_spacing;
    int total_h = num_rows * y_spacing;

    // Build combined grid
    std::vector<std::vector<bool>> combined(total_h, std::vector<bool>(total_w, false));

    for (int i = 0; i < n; i++) {
        int col = i % SUMMARY_COLS;
        int row = i / SUMMARY_COLS;
        int base_x = col * x_spacing;
        int base_y = row * y_spacing;

        const auto& sol = solutions[i];
        for (int y = 0; y < sol.height; y++)
            for (int x = 0; x < sol.width; x++)
                if (sol.cells[y][x])
                    combined[base_y + y][base_x + x] = true;
    }

    // Write RLE header
    out << "x = " << total_w << ", y = " << total_h << ", rule = B3/S23\n";

    // Encode RLE
    int pending_eol = 0;
    for (int y = 0; y < total_h; y++) {
        // Find rightmost alive cell
        int last_alive = -1;
        for (int x = total_w - 1; x >= 0; x--) {
            if (combined[y][x]) { last_alive = x; break; }
        }

        if (last_alive == -1) {
            pending_eol++;
            continue;
        }

        // Flush pending EOLs
        if (pending_eol > 0) {
            if (pending_eol > 1) out << pending_eol;
            out << "$";
            pending_eol = 0;
        }

        // Encode row up to last alive cell
        int x = 0;
        while (x <= last_alive) {
            bool state = combined[y][x];
            int run = 1;
            while (x + run <= last_alive && combined[y][x + run] == state) run++;

            // Skip trailing dead cells (shouldn't happen since we stop at last_alive)
            if (!state && x + run > last_alive) break;

            if (run > 1) out << run;
            out << (state ? 'o' : 'b');
            x += run;
        }

        pending_eol = 1;
    }

    out << "!\n";
}

// Write LifeHistory debug RLE showing relevant cells.
// State 0 = dead non-relevant, 1 = alive non-relevant,
// 2 (B) = dead relevant, 3 (C) = alive relevant.
inline void write_debug_rle(const std::vector<FullSolution>& solutions,
                             std::ostream& out) {
    if (solutions.empty()) return;

    int cell_w = solutions[0].width;
    int cell_h = solutions[0].height;
    int x_spacing = round_up_to(cell_w + SUMMARY_CELL_BUFFER, 10);
    int y_spacing = round_up_to(cell_h + SUMMARY_CELL_BUFFER, 10);

    int n = (int)solutions.size();
    int num_cols = std::min(n, SUMMARY_COLS);
    int num_rows = (n + SUMMARY_COLS - 1) / SUMMARY_COLS;

    int total_w = num_cols * x_spacing;
    int total_h = num_rows * y_spacing;

    // Build combined grid: 0=dead, 1=alive, 2=dead+relevant, 3=alive+relevant
    std::vector<std::vector<int>> combined(total_h, std::vector<int>(total_w, 0));

    for (int i = 0; i < n; i++) {
        int col = i % SUMMARY_COLS;
        int row = i / SUMMARY_COLS;
        int base_x = col * x_spacing;
        int base_y = row * y_spacing;

        const auto& sol = solutions[i];
        for (int y = 0; y < sol.height; y++)
            for (int x = 0; x < sol.width; x++) {
                int state = 0;
                if (sol.cells[y][x] && sol.relevant[y][x]) state = 3; // C
                else if (!sol.cells[y][x] && sol.relevant[y][x]) state = 2; // B
                else if (sol.cells[y][x]) state = 1; // A/o
                combined[base_y + y][base_x + x] = state;
            }
    }

    out << "x = " << total_w << ", y = " << total_h << ", rule = LifeHistory\n";

    int pending_eol = 0;
    for (int y = 0; y < total_h; y++) {
        // Find rightmost non-zero cell
        int last_nonzero = -1;
        for (int x = total_w - 1; x >= 0; x--) {
            if (combined[y][x] != 0) { last_nonzero = x; break; }
        }

        if (last_nonzero == -1) {
            pending_eol++;
            continue;
        }

        if (pending_eol > 0) {
            if (pending_eol > 1) out << pending_eol;
            out << "$";
            pending_eol = 0;
        }

        // LifeHistory RLE chars: . (state 0), A (state 1), B (state 2), C (state 3)
        static const char state_chars[] = {'.', 'A', 'B', 'C'};
        int x = 0;
        while (x <= last_nonzero) {
            int state = combined[y][x];
            int run = 1;
            while (x + run <= last_nonzero && combined[y][x + run] == state) run++;

            if (state == 0 && x + run > last_nonzero) break; // trailing dead

            if (run > 1) out << run;
            out << state_chars[state];
            x += run;
        }

        pending_eol = 1;
    }

    out << "!\n";
}

// Write summary RLE to file derived from input filename.
inline std::string write_summary_rle_file(const std::vector<FullSolution>& solutions,
                                           const std::string& input_filename) {
    // Derive output filename: replace .toml with _results.rle
    std::string out_name = input_filename;
    auto dot = out_name.rfind('.');
    if (dot != std::string::npos)
        out_name = out_name.substr(0, dot);
    out_name += "_results.rle";

    std::ofstream out(out_name);
    if (!out.is_open()) {
        std::cerr << "Warning: could not open " << out_name << " for writing\n";
        return "";
    }

    write_summary_rle(solutions, out);
    return out_name;
}

// Write debug LifeHistory RLE to file derived from input filename.
inline std::string write_debug_rle_file(const std::vector<FullSolution>& solutions,
                                         const std::string& input_filename) {
    std::string out_name = input_filename;
    auto dot = out_name.rfind('.');
    if (dot != std::string::npos)
        out_name = out_name.substr(0, dot);
    out_name += "_debug.rle";

    std::ofstream out(out_name);
    if (!out.is_open()) {
        std::cerr << "Warning: could not open " << out_name << " for writing\n";
        return "";
    }

    write_debug_rle(solutions, out);
    return out_name;
}

// Generate blocking clause: negate this exact catalyst assignment
inline std::vector<int> blocking_clause(const Grid& grid, const SolverResult& result) {
    std::vector<int> clause;
    for (auto [wx, wy] : grid.catalyst_positions) {
        int var = grid.catalyst_var_at(wx, wy);
        if (var < 2) continue;
        int sat_var = var - 1;
        bool alive = result.solution.count(sat_var) > 0;
        // Negate: if alive, add -var; if dead, add +var
        clause.push_back(alive ? -sat_var : sat_var);
    }
    return clause;
}
