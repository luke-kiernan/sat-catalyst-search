#pragma once
#include <string>
#include <vector>
#include <stdexcept>

// LifeHistory cell states. State numbers match Barrister's convention so
// inputs are interchangeable (see Barrister README):
//   1 active pattern, 2 unknown search region,
//   3 catalyst ON cell (initially ON, may be perturbed, must recover),
//   4 catalyst ON cell that is initially OFF (must recover to ON),
//   5 "stator" ON cell that must not be disturbed.
enum class CellState : int {
    DEAD = 0,         // . or b — background dead
    ACTIVE = 1,       // A or o — active pattern (state 1)
    UNKNOWN = 2,      // B — unknown/search region (state 2)
    NON_STATOR = 3,   // C — ON catalyst cell, initially ON, may be active, must recover (state 3)
    INIT_OFF = 4,     // D — catalyst cell ON in stable form but OFF at gen 0, must recover (state 4)
    STATOR = 5,       // E — ON catalyst cell that must not be disturbed (state 5)
};

struct ParsedRLE {
    int width = 0;
    int height = 0;
    std::vector<std::vector<CellState>> grid;  // [y][x]
};

// Parse a LifeHistory RLE string.
// Supports: b/. = dead, o/A = active, B = unknown, C = stator, D = unused, E = non-stator
// Run-length encoding with digits, $ = newline, ! = end
inline ParsedRLE parse_rle(const std::string& rle) {
    ParsedRLE result;

    // Find header line: "x = W, y = H, rule = ..."
    size_t pos = 0;

    // Skip comment lines (starting with #)
    while (pos < rle.size() && rle[pos] == '#') {
        pos = rle.find('\n', pos);
        if (pos == std::string::npos) break;
        pos++;
    }

    // Skip whitespace
    while (pos < rle.size() && (rle[pos] == ' ' || rle[pos] == '\n' || rle[pos] == '\r'))
        pos++;

    // Parse header
    if (pos < rle.size() && rle[pos] == 'x') {
        // Parse "x = W, y = H, rule = ..."
        auto skip_ws = [&]() {
            while (pos < rle.size() && (rle[pos] == ' ' || rle[pos] == '\t')) pos++;
        };
        auto parse_int = [&]() -> int {
            int val = 0;
            while (pos < rle.size() && rle[pos] >= '0' && rle[pos] <= '9') {
                val = val * 10 + (rle[pos] - '0');
                pos++;
            }
            return val;
        };

        // x = W
        pos++; skip_ws();
        if (pos < rle.size() && rle[pos] == '=') pos++;
        skip_ws();
        result.width = parse_int();

        // , y = H
        skip_ws();
        if (pos < rle.size() && rle[pos] == ',') pos++;
        skip_ws();
        if (pos < rle.size() && rle[pos] == 'y') pos++;
        skip_ws();
        if (pos < rle.size() && rle[pos] == '=') pos++;
        skip_ws();
        result.height = parse_int();

        // Skip rest of header line (rule = ...)
        while (pos < rle.size() && rle[pos] != '\n') pos++;
        if (pos < rle.size()) pos++;
    }

    if (result.width <= 0 || result.height <= 0) {
        throw std::runtime_error("RLE: invalid or missing header (width=" +
            std::to_string(result.width) + ", height=" + std::to_string(result.height) + ")");
    }

    // Initialize grid
    result.grid.resize(result.height, std::vector<CellState>(result.width, CellState::DEAD));

    // Parse pattern data
    int x = 0, y = 0;
    int run_count = 0;

    while (pos < rle.size()) {
        char c = rle[pos++];

        if (c >= '0' && c <= '9') {
            run_count = run_count * 10 + (c - '0');
            continue;
        }

        int count = (run_count == 0) ? 1 : run_count;
        run_count = 0;

        CellState state;
        bool is_cell = true;

        switch (c) {
            case 'b': case '.':
                state = CellState::DEAD;
                break;
            case 'o': case 'A':
                state = CellState::ACTIVE;
                break;
            case 'B':
                state = CellState::UNKNOWN;
                break;
            case 'C':
                state = CellState::NON_STATOR;
                break;
            case 'D':
                state = CellState::INIT_OFF;
                break;
            case 'E':
                state = CellState::STATOR;
                break;
            case '$':
                y += count;
                x = 0;
                is_cell = false;
                break;
            case '!':
                return result;
            case '\n': case '\r': case ' ':
                is_cell = false;
                break;
            default:
                is_cell = false;
                break;
        }

        if (is_cell) {
            for (int i = 0; i < count; i++) {
                if (y < result.height && x < result.width) {
                    result.grid[y][x] = state;
                }
                x++;
            }
        }
    }

    return result;
}
