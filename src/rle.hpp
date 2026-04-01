#pragma once
#include <string>
#include <vector>
#include <stdexcept>

// LifeHistory cell states
enum class CellState : int {
    DEAD = 0,       // . or b — background dead
    ACTIVE = 1,     // A or o — active pattern (state 1)
    UNKNOWN = 2,    // B — unknown/search region (state 2)
    FORCED_ON = 3,  // C — forced ON (state 3)
    FORCED_ON_INIT_OFF = 4, // D — forced ON initially OFF (state 4)
    // State 5 (E) = stator, deferred to v2
};

struct ParsedRLE {
    int width = 0;
    int height = 0;
    std::vector<std::vector<CellState>> grid;  // [y][x]
};

// Parse a LifeHistory RLE string.
// Supports: b/. = dead, o/A = active, B = unknown, C = forced ON, D = forced ON init OFF
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
                state = CellState::FORCED_ON;
                break;
            case 'D':
                state = CellState::FORCED_ON_INIT_OFF;
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
