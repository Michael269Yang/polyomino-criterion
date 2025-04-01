#ifndef ISOHEDRAL_H
#define ISOHEDRAL_H

#include "boundary.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <map>
#include <utility>
#include <vector>

// TODO: Put this stuff in a namespace.

using Factor = std::pair<int, int>;

std::vector<Factor> is_double_palindrome(const Factor& F, const std::vector<std::vector<Factor>>& palindrome_factor_starts, const std::vector<std::vector<Factor>>& palindrome_factor_ends, int n);


struct IsohedralChecker {

int minAngle = 90;

std::map<std::pair<int, int>, std::pair<int, int>> COMPLEMENT = {
  {{0, 1}, {0, -1}},
  {{0, -1}, {0, 1}},
  {{1, 0}, {-1, 0}},
  {{-1, 0}, {1, 0}}
};

std::map<std::pair<int, int>, std::pair<int, int>> CCW = {
  {{0, 1}, {-1, 0}},
  {{1, 0}, {0, 1}},
  {{0, -1}, {1, 0}},
  {{-1, 0}, {0, -1}}
};

std::map<std::pair<int, int>, std::pair<int, int>> CW = {
    {{-1, 0}, {0, 1}},
    {{0, 1}, {1, 0}},
    {{1, 0}, {0, -1}},
    {{0, -1}, {-1, 0}}
};

std::map<int, std::map<std::pair<int, int>, std::pair<int, int>>> REFL = {
  {-45, {{{0, 1}, {-1, 0}}, {{1, 0}, {0, -1}}, {{0, -1}, {1, 0}}, {{-1, 0}, {0, 1}}}},
  {0, {{{0, 1}, {0, -1}}, {{1, 0}, {1, 0}}, {{0, -1}, {0, 1}}, {{-1, 0}, {-1, 0}}}},
  {45, {{{0, 1}, {1, 0}}, {{1, 0}, {0, 1}}, {{0, -1}, {-1, 0}}, {{-1, 0}, {0, -1}}}},
  {90, {{{0, 1}, {0, 1}}, {{1, 0}, {-1, 0}}, {{0, -1}, {0, -1}}, {{-1, 0}, {1, 0}}}},
};

std::pair<int, int> iteratedCcw(std::pair<int, int> dir, int numIters);

bool is_reflect_square_factor(const boundaryword& P, int i, int j, int theta);

boundaryword inv_comp(const boundaryword& S);

std::vector<Factor> admissible_mirror_factors(const boundaryword& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_mirror_factor_pairs(const boundaryword& P);

std::vector<Factor> admissible_rotation_factors(const boundaryword& P, int theta);

std::vector<Factor> admissible_rotadrome_factors(const boundaryword& P, int theta);

std::vector<Factor> admissible_reflect_square_factors(const boundaryword& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_reflect_square_factor_pairs(const boundaryword& P, int theta);

std::vector<Factor> has_translation_tiling(const boundaryword& P);

std::vector<Factor> has_half_turn_tiling(const boundaryword& P);

std::vector<Factor> has_quarter_turn_tiling(const boundaryword& P);

std::vector<Factor> has_type_1_reflection_tiling(const boundaryword& P);

std::vector<Factor> has_type_2_reflection_tiling(const boundaryword& P);

std::vector<Factor> has_type_1_half_turn_reflection_tiling(const boundaryword& P);

std::vector<Factor> has_type_2_half_turn_reflection_tiling(const boundaryword& P);

std::vector<Factor> has_case_7_tiling(const boundaryword& P);

std::vector<Factor> has_case_8a_tiling(const boundaryword& P);

std::vector<Factor> has_case_8b_tiling(const boundaryword& P);


bool has_isohedral_tiling(const boundaryword& P);

};

#endif // ISOHEDRAL_H

