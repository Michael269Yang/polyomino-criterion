#ifndef ISOHEDRAL_H
#define ISOHEDRAL_H

#include "boundary.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <map>
#include <utility>
#include <vector>

constexpr int MAX_BND = 50;

using Factor = std::pair<int, int>;

// A wrapper around stl array that provides range based iteration
// over only the filled elements.
template <typename T, size_t N>
struct PartialArray {
  std::array<T, N> data;
  size_t filled_count;

  T& operator[](size_t index) {
    return data[index];
  }

  auto begin() const { 
    return data.cbegin();
  }

  auto end() const {
    return data.cbegin() + filled_count;
  }
};

template <size_t N>
using FactorArray = PartialArray<Factor, N>;

using FactorPair = std::pair<Factor, Factor>;

template <size_t N>
using FactorPairArray = PartialArray<FactorPair, N>;

// Returns whether factor F = [start, end] (both inclusive) is a 
// double palidrome. Assumes palindrome_factor_starts, palindrome_factor_ends
// are precomputed with the start and ends of the maximal palindrome factors
// in the boundary word.
bool is_double_palindrome(const Factor& F, const std::vector<std::vector<Factor>>& palindrome_factor_starts, const std::vector<std::vector<Factor>>& palindrome_factor_ends, int n);

// A struct which checks whether a boundary word corresponds to 
// a shape that tiles the plane isohedrally.
//
// The default fields of the struct are for the polyomino grid.
// Assumes equal angle between edges in the grid.
struct IsohedralChecker {

// The angle between edges of the grid.
int minAngle = 90;

// The complement of edges in the grid.
std::map<std::pair<int, int>, std::pair<int, int>> COMPLEMENT = {
  {{0, 1}, {0, -1}},
  {{0, -1}, {0, 1}},
  {{1, 0}, {-1, 0}},
  {{-1, 0}, {1, 0}}
};

// The result of rotating an edge counterclockwise.
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

// The results of reflecting edges in the grid.
std::map<int, std::map<std::pair<int, int>, std::pair<int, int>>> REFL = {
  {-45, {{{0, 1}, {-1, 0}}, {{1, 0}, {0, -1}}, {{0, -1}, {1, 0}}, {{-1, 0}, {0, 1}}}},
  {0, {{{0, 1}, {0, -1}}, {{1, 0}, {1, 0}}, {{0, -1}, {0, 1}}, {{-1, 0}, {-1, 0}}}},
  {45, {{{0, 1}, {1, 0}}, {{1, 0}, {0, 1}}, {{0, -1}, {-1, 0}}, {{-1, 0}, {0, -1}}}},
  {90, {{{0, 1}, {0, 1}}, {{1, 0}, {-1, 0}}, {{0, -1}, {0, -1}}, {{-1, 0}, {1, 0}}}},
};

// Rotates dir counterclockwise numIters times.
std::pair<int, int> iteratedCcw(std::pair<int, int> dir, int numIters);

// Returns if P[i, j] is of the form A refl(A).
bool is_reflect_square_factor(const boundaryword& P, int i, int j, int theta);

// Reverses S and takes the componentwise complement.
boundaryword inv_comp(const boundaryword& S);

FactorArray<2*MAX_BND> admissible_mirror_factors(const boundaryword& P);

FactorPairArray<MAX_BND*MAX_BND> admissible_gapped_mirror_factor_pairs(const boundaryword& P);

FactorArray<2*MAX_BND> admissible_rotadrome_factors(const boundaryword& P, int theta);

FactorArray<MAX_BND*MAX_BND> admissible_reflect_square_factors(const boundaryword& P);

std::vector<std::pair<Factor, Factor>> admissible_gapped_reflect_square_factor_pairs(const boundaryword& P, int theta);

bool has_translation_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& mirror_factors);

bool has_half_turn_tiling(const boundaryword& P, const FactorPairArray<MAX_BND*MAX_BND>& mirror_factor_pairs, const FactorArray<2*MAX_BND>& palin_factors);

bool has_quarter_turn_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& ninety_factors, const FactorArray<2*MAX_BND>& palin_factors);

bool has_type_1_reflection_tiling(const boundaryword& P, const FactorArray<MAX_BND*MAX_BND>& reflect_square_factors, const FactorPairArray<MAX_BND*MAX_BND>& mirror_factor_pairs);

bool has_type_2_reflection_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& mirror_factors);

bool has_type_1_half_turn_reflection_tiling(const boundaryword& P, const FactorPairArray<MAX_BND*MAX_BND>& partial_mirror_factor_pairs, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<MAX_BND*MAX_BND>& reflect_square_factors);

bool has_type_2_half_turn_reflection_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors);

bool has_case_7_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& onetwenty_factors);

bool has_case_8a_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<2*MAX_BND>& sixty_factors, const FactorArray<2*MAX_BND>& onetwenty_factors);

bool has_case_8b_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<2*MAX_BND>& sixty_factors, const FactorArray<2*MAX_BND>& onetwenty_factors);

bool has_isohedral_tiling(const boundaryword& P);

};

#endif // ISOHEDRAL_H

