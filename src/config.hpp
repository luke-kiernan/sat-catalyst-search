#pragma once
#include <string>
#include <array>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include "toml.hpp"
#include "rle.hpp"
#include "rule.hpp"

struct SearchConfig {
    // Parsed pattern
    ParsedRLE pattern;
    int center_x = 0;
    int center_y = 0;

    // Cellular automaton rule (default B3/S23 for backward compatibility).
    Rule rule = parse_rule("B3/S23");

    // Timing
    std::array<int, 2> first_active_range = {0, 0};
    std::array<int, 2> active_window_range = {0, 0};
    int min_stable_interval = 0;

    // Perturbation limits
    int max_active_cells = -1;       // -1 = no limit
    int max_ever_active_cells = -1;  // -1 = no limit
    std::array<int, 2> ever_active_bounds = {-1, -1};  // -1 = no limit

    // Computed
    int total_generations() const {
        return first_active_range[1] + active_window_range[1] + min_stable_interval;
    }
};

inline SearchConfig parse_config(const std::string& filename) {
    SearchConfig config;

    auto data = toml::parse(filename);

    // Pattern
    std::string rle_str = toml::find<std::string>(data, "pattern");
    config.pattern = parse_rle(rle_str);

    // Pattern center
    auto center = toml::find<std::vector<int>>(data, "pattern-center");
    if (center.size() != 2)
        throw std::runtime_error("pattern-center must be [x, y]");
    config.center_x = center[0];
    config.center_y = center[1];

    // Timing parameters
    auto far = toml::find<std::vector<int>>(data, "first-active-range");
    if (far.size() != 2)
        throw std::runtime_error("first-active-range must be [min, max]");
    config.first_active_range = {far[0], far[1]};

    auto awr = toml::find<std::vector<int>>(data, "active-window-range");
    if (awr.size() != 2)
        throw std::runtime_error("active-window-range must be [min, max]");
    config.active_window_range = {awr[0], awr[1]};

    config.min_stable_interval = toml::find<int>(data, "min-stable-interval");

    // Optional rule (defaults to B3/S23).
    if (data.contains("rule"))
        config.rule = parse_rule(toml::find<std::string>(data, "rule"));

    // Optional perturbation limits
    if (data.contains("max-active-cells"))
        config.max_active_cells = toml::find<int>(data, "max-active-cells");

    if (data.contains("max-ever-active-cells"))
        config.max_ever_active_cells = toml::find<int>(data, "max-ever-active-cells");

    if (data.contains("ever-active-bounds")) {
        auto eab = toml::find<std::vector<int>>(data, "ever-active-bounds");
        if (eab.size() != 2)
            throw std::runtime_error("ever-active-bounds must be [w, h]");
        config.ever_active_bounds = {eab[0], eab[1]};
    }

    return config;
}
