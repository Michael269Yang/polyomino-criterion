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

bool is_reflect_square_factor(const std::string& P, int i, int j, int theta);

std::string inv_comp(const std::string& S);

std::vector<Factor> admissible_mirror_factors(const std::string& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_mirror_factor_pairs(const std::string& P);

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

bool has_isohedral_tiling(const std::string& P);

#endif // ISOHEDRAL_H

