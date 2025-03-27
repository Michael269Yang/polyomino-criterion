#ifndef ISOHEDRAL_H
#define ISOHEDRAL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// TODO: Put this stuff in a namespace.

using Factor = std::pair<int, int>;

std::vector<Factor> is_double_palindrome(const Factor& F, const std::vector<std::vector<Factor>>& palindrome_factor_starts, const std::vector<std::vector<Factor>>& palindrome_factor_ends, int n);


struct IsohedralChecker {

int minAngle = 90;

std::unordered_map<char, char> COMPLEMENT = {
  {'N', 'S'},
  {'S', 'N'},
  {'E', 'W'},
  {'W', 'E'}
};

std::unordered_map<char, char> CCW = {
  {'N', 'W'},
  {'E', 'N'},
  {'S', 'E'},
  {'W', 'S'}
};

std::unordered_map<char, char> CW = {
    {'W', 'N'},
    {'N', 'E'},
    {'E', 'S'},
    {'S', 'W'}
};

std::unordered_map<int, std::unordered_map<char, char>> REFL = {
  {-45, {{'N', 'W'}, {'E', 'S'}, {'S', 'E'}, {'W', 'N'}}},
  {0, {{'N', 'S'}, {'E', 'E'}, {'S', 'N'}, {'W', 'W'}}},
  {45, {{'N', 'E'}, {'E', 'N'}, {'S', 'W'}, {'W', 'S'}}},
  {90, {{'N', 'N'}, {'E', 'W'}, {'S', 'S'}, {'W', 'E'}}},
};

char iteratedCw(char dir, int numIters);

bool is_reflect_square_factor(const std::string& P, int i, int j, int theta);

std::string inv_comp(const std::string& S);

std::vector<Factor> admissible_mirror_factors(const std::string& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_mirror_factor_pairs(const std::string& P);

std::vector<Factor> admissible_rotation_factors(const std::string& P, int theta);

std::vector<Factor> admissible_rotadrome_factors(const std::string& P, int theta);

std::vector<Factor> admissible_reflect_square_factors(const std::string& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_reflect_square_factor_pairs(const std::string& P, int theta);

std::vector<Factor> has_translation_tiling(const std::string& P);

std::vector<Factor> has_half_turn_tiling(const std::string& P);

std::vector<Factor> has_quarter_turn_tiling(const std::string& P);

std::vector<Factor> has_type_1_reflection_tiling(const std::string& P);

std::vector<Factor> has_type_2_reflection_tiling(const std::string& P);

std::vector<Factor> has_type_1_half_turn_reflection_tiling(const std::string& P);

std::vector<Factor> has_type_2_half_turn_reflection_tiling(const std::string& P);

std::vector<Factor> has_case_7_tiling(const std::string& P);

std::vector<Factor> has_case_8a_tiling(const std::string& P);

std::vector<Factor> has_case_8b_tiling(const std::string& P);


bool has_isohedral_tiling(const std::string& P);

};

#endif // ISOHEDRAL_H

