#pragma once
#include <array>
#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Isotropic Hensel-notation rule (2-state, 3x3 Moore neighborhood).
//
// 9-bit neighborhood encoding (matches hensel.py's get_9bit):
//   bit 0 1 2     NW  N  NE
//   bit 3 4 5     W   C  E
//   bit 6 7 8     SW  S  SE
//
// Rule::alive_next[n] = true iff a cell with 9-bit neighborhood n
// (including its own current state at bit 4) becomes alive at t+1.
struct Rule {
    std::array<bool, 512> alive_next{};
    std::string canonical;  // e.g. "B3/S23" — preserved for RLE output

    bool evolves_to(int nine_bit) const { return alive_next[nine_bit]; }

    // A still-life configuration: the rule's next state equals the center bit.
    bool is_stable(int nine_bit) const {
        bool center = ((nine_bit >> 4) & 1) != 0;
        return alive_next[nine_bit] == center;
    }
};

namespace hensel_detail {

// Letters valid for each digit (matches hensel.py's conditiondict).
// Empty string means "the bare digit, no letter suffix".
inline const std::map<char, std::vector<std::string>>& letter_table() {
    static const std::map<char, std::vector<std::string>> t = {
        {'0', {""}},
        {'1', {"c","e"}},
        {'2', {"a","c","e","i","k","n"}},
        {'3', {"a","c","e","i","j","k","n","q","r","y"}},
        {'4', {"a","c","e","i","j","k","n","q","r","t","w","y","z"}},
        {'5', {"a","c","e","i","j","k","n","q","r","y"}},
        {'6', {"a","c","e","i","k","n"}},
        {'7', {"c","e"}},
        {'8', {""}},
    };
    return t;
}

// Direct port of hensel.py's get_9bit: maps a "B<digit><letter>" condition
// to the list of 9-bit neighborhood values (center bit always 0 here;
// survival conditions are derived by adding 16 = bit 4).
inline const std::map<std::string, std::vector<int>>& nine_bit_table() {
    static const std::map<std::string, std::vector<int>> t = {
        {"B0",  {0}},
        {"B1c", {1, 4, 64, 256}},
        {"B1e", {2, 8, 32, 128}},
        {"B2a", {3, 6, 9, 36, 72, 192, 288, 384}},
        {"B2c", {5, 65, 260, 320}},
        {"B2e", {10, 34, 136, 160}},
        {"B2i", {40, 130}},
        {"B2k", {12, 33, 66, 96, 129, 132, 258, 264}},
        {"B2n", {68, 257}},
        {"B3a", {11, 38, 200, 416}},
        {"B3c", {69, 261, 321, 324}},
        {"B3e", {42, 138, 162, 168}},
        {"B3i", {7, 73, 292, 448}},
        {"B3j", {14, 35, 74, 137, 164, 224, 290, 392}},
        {"B3k", {98, 140, 161, 266}},
        {"B3n", {13, 37, 67, 193, 262, 328, 352, 388}},
        {"B3q", {70, 76, 100, 196, 259, 265, 289, 385}},
        {"B3r", {41, 44, 104, 131, 134, 194, 296, 386}},
        {"B3y", {97, 133, 268, 322}},
        {"B4a", {15, 39, 75, 201, 294, 420, 456, 480}},
        {"B4c", {325}},
        {"B4e", {170}},
        {"B4i", {45, 195, 360, 390}},
        {"B4j", {106, 142, 163, 169, 172, 226, 298, 394}},
        {"B4k", {99, 141, 165, 225, 270, 330, 354, 396}},
        {"B4n", {71, 77, 263, 293, 329, 356, 449, 452}},
        {"B4q", {102, 204, 267, 417}},
        {"B4r", {43, 46, 139, 166, 202, 232, 418, 424}},
        {"B4t", {105, 135, 300, 450}},
        {"B4w", {78, 228, 291, 393}},
        {"B4y", {101, 197, 269, 323, 326, 332, 353, 389}},
        {"B4z", {108, 198, 297, 387}},
        {"B5a", {79, 295, 457, 484}},
        {"B5c", {171, 174, 234, 426}},
        {"B5e", {327, 333, 357, 453}},
        {"B5i", {47, 203, 422, 488}},
        {"B5j", {103, 205, 271, 331, 358, 421, 460, 481}},
        {"B5k", {229, 334, 355, 397}},
        {"B5n", {107, 143, 167, 233, 302, 428, 458, 482}},
        {"B5q", {110, 206, 230, 236, 299, 395, 419, 425}},
        {"B5r", {109, 199, 301, 361, 364, 391, 451, 454}},
        {"B5y", {173, 227, 362, 398}},
        {"B6a", {111, 207, 303, 423, 459, 486, 489, 492}},
        {"B6c", {175, 235, 430, 490}},
        {"B6e", {335, 359, 461, 485}},
        {"B6i", {365, 455}},
        {"B6k", {231, 237, 363, 366, 399, 429, 462, 483}},
        {"B6n", {238, 427}},
        {"B7c", {239, 431, 491, 494}},
        {"B7e", {367, 463, 487, 493}},
        {"B8",  {495}},
    };
    return t;
}

// Resolve a "B<digit><letters?>" or "S<...>" condition into its 9-bit values.
// For S conditions, defers to B with +16 (sets the center bit) — matches hensel.py.
inline std::vector<int> nine_bit_for(const std::string& cond) {
    if (cond.empty()) return {};
    if (cond[0] == 'S') {
        std::string b = "B" + cond.substr(1);
        auto vs = nine_bit_for(b);
        for (auto& v : vs) v += 16;
        return vs;
    }
    auto& tbl = nine_bit_table();
    auto it = tbl.find(cond);
    if (it == tbl.end())
        throw std::runtime_error("Rule: unknown condition '" + cond + "'");
    return it->second;
}

// Port of hensel.py's parserule. Returns sub-conditions like
// {"B3a","B3c",...,"S2a",...}. Handles bare digits ("B3" → all 10 letters),
// explicit letter lists ("B3aei"), and negation ("B3-jr" → all except j,r).
inline std::vector<std::string> parse_conditions(std::string rule) {
    std::transform(rule.begin(), rule.end(), rule.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    rule.erase(std::remove(rule.begin(), rule.end(), '/'), rule.end());

    if (rule.empty() || rule[0] != 'b')
        throw std::runtime_error("Rule: must start with 'B' or 'b'");
    if (rule.find('s') == std::string::npos)
        throw std::runtime_error("Rule: missing 'S' section");

    std::vector<std::string> conditions;
    bool birth = true;
    std::string cstring;
    rule.push_back('|');  // sentinel to flush the final segment

    auto flush = [&](char prefix_letter) {
        // cstring at this point is like "3", "3a", "3aei", "3-jr"
        if (cstring.empty()) return;
        char digit = cstring[0];
        if (digit < '0' || digit > '8')
            throw std::runtime_error("Rule: bad digit in '" + cstring + "'");
        char tag = static_cast<char>(std::toupper(prefix_letter));
        const auto& letters = letter_table().at(digit);

        if (cstring.size() == 1) {
            // Bare digit → all letters for that digit.
            for (const auto& L : letters)
                conditions.push_back(std::string(1, tag) + digit + L);
        } else if (cstring[1] == '-') {
            // Negation: all letters EXCEPT those listed after '-'.
            std::vector<std::string> remaining = letters;
            for (size_t i = 2; i < cstring.size(); i++) {
                std::string lit(1, cstring[i]);
                remaining.erase(std::remove(remaining.begin(), remaining.end(), lit),
                                remaining.end());
            }
            for (const auto& L : remaining)
                conditions.push_back(std::string(1, tag) + digit + L);
        } else {
            // Explicit letter list: keep only listed letters that are valid.
            for (size_t i = 1; i < cstring.size(); i++) {
                std::string lit(1, cstring[i]);
                if (std::find(letters.begin(), letters.end(), lit) != letters.end())
                    conditions.push_back(std::string(1, tag) + digit + lit);
            }
        }
    };

    char current_tag = 'b';  // 'b' until we hit 's'
    for (char ch : rule) {
        bool is_digit = (ch >= '0' && ch <= '8');
        if (is_digit || ch == 'b' || ch == 's' || ch == '|') {
            // Boundary — flush accumulated cstring under whatever section we were in.
            flush(current_tag);
            cstring.clear();
            if (ch == 'b') { birth = true; current_tag = 'b'; continue; }
            if (ch == 's') { birth = false; current_tag = 's'; continue; }
            if (ch == '|') continue;
            // digit: start a new cstring
            cstring.push_back(ch);
        } else {
            // letter or '-' — accumulate
            cstring.push_back(ch);
        }
    }
    (void)birth;
    return conditions;
}

}  // namespace hensel_detail

inline Rule parse_rule(const std::string& s) {
    Rule r;
    auto conds = hensel_detail::parse_conditions(s);
    for (const auto& c : conds)
        for (int v : hensel_detail::nine_bit_for(c))
            r.alive_next[v] = true;
    r.canonical = s;
    return r;
}
