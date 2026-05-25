#include <iostream>
#include <cassert>
#include <string>
#include "rule.hpp"

// Pack a 3x3 neighborhood into a 9-bit int with the layout:
//   bit 0 1 2     NW  N  NE
//   bit 3 4 5     W   C  E
//   bit 6 7 8     SW  S  SE
static int pack(bool nw, bool n, bool ne, bool w, bool c, bool e, bool sw, bool s, bool se) {
    return nw | (n<<1) | (ne<<2) | (w<<3) | (c<<4) | (e<<5) | (sw<<6) | (s<<7) | (se<<8);
}

// 1. Round-trip the canonical CGoL rule and verify the well-known transitions.
//    Counts cross-check: B3/S23 = C(8,3) + C(8,2) + C(8,3) = 56 + 28 + 56 = 140.
static void test_b3s23() {
    std::cout << "-- B3/S23 transitions and count\n";
    Rule r = parse_rule("B3/S23");

    int count = 0;
    for (int n = 0; n < 512; n++) if (r.alive_next[n]) count++;
    assert(count == 140);

    // Dead cell with 3 alive neighbors → birth
    assert(r.evolves_to(pack(1,1,1, 0,0,0, 0,0,0)) == true);
    // Dead cell with 2 alive neighbors → no birth
    assert(r.evolves_to(pack(1,1,0, 0,0,0, 0,0,0)) == false);
    // Alive cell with 2 alive neighbors → survive
    assert(r.evolves_to(pack(1,1,0, 0,1,0, 0,0,0)) == true);
    // Alive cell with 1 alive neighbor → die
    assert(r.evolves_to(pack(1,0,0, 0,1,0, 0,0,0)) == false);
    // Alive cell with 4 alive neighbors → die
    assert(r.evolves_to(pack(1,1,1, 1,1,0, 0,0,0)) == false);
}

// 2. Stability: block (2x2) is a still life, blinker is not.
static void test_stability() {
    std::cout << "-- is_stable distinguishes still lifes from oscillators\n";
    Rule r = parse_rule("B3/S23");

    // Block — center cell of a 2x2 with 3 alive neighbors (the other 3 corners).
    // Center=1, neighbors=3 → next state alive → matches center → stable.
    assert(r.is_stable(pack(1,1,0, 1,1,0, 0,0,0)) == true);

    // Blinker horizontal: center alive with 2 alive neighbors (left+right).
    // Next state alive → stable as center, BUT the cell above the blinker
    // (dead with 3 alive neighbors) is NOT stable — it births.
    assert(r.is_stable(pack(0,0,0, 1,1,1, 0,0,0)) == true);
    assert(r.is_stable(pack(1,1,1, 0,0,0, 0,0,0)) == false);  // dead, 3 nbrs → births
}

// 3. Hensel sub-conditions: explicit letter list narrows which patterns fire.
//    Use the rule's notation to drop one specific 3-neighbor pattern shape
//    and verify the right 9-bit values are/aren't in alive_next.
static void test_isotropic_letter_list() {
    std::cout << "-- Hensel letter lists restrict transitions correctly\n";
    // B3a only: just the 'a' patterns of B3 → 4 of the 56 total B3 patterns.
    // hensel.py: B3a = [11, 38, 200, 416].
    Rule r = parse_rule("B3a/S23");

    int b_count = 0;
    for (int n = 0; n < 512; n++)
        if (r.alive_next[n] && ((n >> 4) & 1) == 0) b_count++;
    assert(b_count == 4);

    assert(r.evolves_to(11) == true);
    // A 3-neighbor pattern that is NOT B3a — e.g. B3i = [7, 73, 292, 448].
    assert(r.evolves_to(7) == false);
}

// 4. Negation syntax: "B3-jr" means all B3 letters EXCEPT j and r.
//    Total: 10 - 2 = 8 letters → 56 - 8 - 8 = 40 patterns (B3j=8, B3r=8).
static void test_negation() {
    std::cout << "-- 'B<digit>-<letters>' syntax excludes the listed letters\n";
    Rule r = parse_rule("B3-jr/S23");

    int b3_count = 0;
    for (int n = 0; n < 512; n++)
        if (r.alive_next[n] && ((n >> 4) & 1) == 0) b3_count++;
    assert(b3_count == 56 - 8 - 8);

    // B3j sample: 14 should NOT fire; B3a sample: 11 should fire.
    assert(r.evolves_to(14) == false);
    assert(r.evolves_to(11) == true);
}

// 5. Canonical string is preserved (used for RLE output later).
static void test_canonical_preserved() {
    std::cout << "-- canonical string round-trips\n";
    Rule r = parse_rule("B36/S23");
    assert(r.canonical == "B36/S23");
}

int main() {
    std::cout << "=== Rule (Hensel) tests ===\n";
    test_b3s23();
    test_stability();
    test_isotropic_letter_list();
    test_negation();
    test_canonical_preserved();
    std::cout << "=== All rule tests passed ===\n";
    return 0;
}
