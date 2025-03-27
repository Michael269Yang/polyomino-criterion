#include "isohedral.h"

#include <string>

using namespace std;

/****************************************************
* Code for tests 
* This section contains:
*   - Helper functions 
*   - Functions that test functionality.
*   - A main function which runs our test suit. 
****************************************************/

IsohedralChecker ominoChecker;

string create_square(int i) {
  string square(4 * i, ' ');
  fill(square.begin(), square.begin() + i, 'N');
  fill(square.begin() + i, square.begin() + 2 * i, 'E');
  fill(square.begin() + 2 * i, square.begin() + 3 * i, 'S');
  fill(square.begin() + 3 * i, square.end(), 'W');
  return square;
}

string create_long_rectangle(int i) {
  string rectangle(2 * i + 2, ' ');
  rectangle[0] = 'N';
  fill(rectangle.begin() + 1, rectangle.begin() + i + 1, 'E');
  rectangle[i + 1] = 'S';
  fill(rectangle.begin() + i + 2, rectangle.end(), 'W');
  return rectangle;
}

string create_tall_rectangle(int i) {
  string rectangle(2 * i + 2, ' ');
  fill(rectangle.begin(), rectangle.begin() + i, 'N');
  rectangle[i] = 'E';
  fill(rectangle.begin() + i + 1, rectangle.begin() + 2 * i + 1, 'S');
  rectangle[2 * i + 1] = 'W';
  return rectangle;
}

string create_2_by_i_rectangle(int i) {
  string rectangle(2 * i + 4, ' ');
  rectangle[0] = 'N';
  rectangle[1] = 'N';
  fill(rectangle.begin() + 2, rectangle.begin() + i + 2, 'E');
  rectangle[i + 2] = 'S';
  rectangle[i + 3] = 'S';
  fill(rectangle.begin() + i + 4, rectangle.end(), 'W');
  return rectangle;
}

string create_rectangle(int n, int m) {
  string rectangle(2 * n + 2 * m, ' ');
  fill(rectangle.begin(), rectangle.begin() + n, 'N');
  fill(rectangle.begin() + n, rectangle.begin() + n + m, 'E');
  fill(rectangle.begin() + n + m, rectangle.begin() + 2*n + m, 'S');
  fill(rectangle.begin() + 2*n + m, rectangle.end(), 'W');
  return rectangle;
}

void test__admissible_rotadrome_factors() {
  vector<Factor> factors = ominoChecker.admissible_rotadrome_factors("NESW", 180);
  set<Factor> factorSet(factors.begin(), factors.end());
  set<Factor> wantSet = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
  assert(factorSet == wantSet);
}

void test__admissible_reflect_square_factors() {
  vector<Factor> factors = ominoChecker.admissible_reflect_square_factors("NEESWW");
  set<Factor> factorSet(factors.begin(), factors.end());
  set<Factor> wantSet = {{0, 1}, {2, 3}, {3, 4}, {5, 0}, {1, 2}, {4, 5}};
  assert(factorSet == wantSet);

  factors = ominoChecker.admissible_reflect_square_factors("NNEESSWW");
  factorSet = set<Factor>(factors.begin(), factors.end());
  wantSet = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 0}, {0, 3}, {2, 5}, {4, 7}, {6, 1}};
  assert(factorSet == wantSet);
}

void test__admissible_gapped_reflect_square_factor_pairs() {
  vector<pair<Factor, Factor>> factors = ominoChecker.admissible_gapped_reflect_square_factor_pairs("NESW", -45);
  set<pair<Factor, Factor>> factorSet(factors.begin(), factors.end());
  set<pair<Factor, Factor>> wantSet = {{{0,0}, {3,3}}, {{1,1}, {2,2}}};
  assert(wantSet == factorSet);

  factors = ominoChecker.admissible_gapped_reflect_square_factor_pairs("NESW", 45);
  factorSet = set<pair<Factor, Factor>>(factors.begin(), factors.end());
  wantSet = {{{0, 0}, {1,1}}, {{2,2}, {3,3}}};
  assert(wantSet == factorSet);

  factors = ominoChecker.admissible_gapped_reflect_square_factor_pairs("NESW", 0);
  factorSet = set<pair<Factor, Factor>>(factors.begin(), factors.end());
  wantSet = {{{0, 0}, {2,2}}};
  assert(wantSet == factorSet);

  factors = ominoChecker.admissible_gapped_reflect_square_factor_pairs("NESW", 90);
  factorSet = set<pair<Factor, Factor>>(factors.begin(), factors.end());
  wantSet = {{{1, 1}, {3,3}}};
  assert(wantSet == factorSet);
}

void test__has_translation_tiling() {
  // Test has_translation_tiling
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_translation_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_translation_tiling(create_long_rectangle(i)).empty());

    // Rectangle in other direction 
    assert(!ominoChecker.has_translation_tiling(create_tall_rectangle(i)).empty());
  }

  assert(!ominoChecker.has_translation_tiling("NNESESSWNW").empty());

  assert(!ominoChecker.has_translation_tiling("NNNESESSWW").empty());
  assert(!ominoChecker.has_translation_tiling("NENENESEESWSWNWSWW").empty());
  assert(ominoChecker.has_translation_tiling("NENENESEESSWWNWSWW").empty());

  assert(ominoChecker.has_translation_tiling("NENNESSSEESWWWNW").empty());
  assert(ominoChecker.has_translation_tiling("WWWNNESENESS").empty());
}

void test__has_half_turn_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_half_turn_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_half_turn_tiling(create_long_rectangle(i)).empty());
    assert(!ominoChecker.has_half_turn_tiling(create_tall_rectangle(i)).empty());
  }

  assert(!ominoChecker.has_half_turn_tiling("NNESESSWWW").empty());
  assert(!ominoChecker.has_half_turn_tiling("NNNESESSWW").empty());

  assert(!ominoChecker.has_half_turn_tiling("NWNEENWNENESESESSWSWNWSW").empty());

  assert(!ominoChecker.has_half_turn_tiling("NENENESEESWSWNWSWW").empty());
  assert(!ominoChecker.has_half_turn_tiling("WWWNNESENESS").empty());

  string B = "ENESEENES";
  string reverseB = B;
  reverse(reverseB.begin(), reverseB.end());
  B += reverseB;

  string C = "SSESSW";
  string reverseC = C;
  reverse(reverseC.begin(), reverseC.end());
  C += reverseC;

  string D = "WWWWWWWWWW";
  
  string E = "NWWNEE";
  string reverseE = E;
  reverse(reverseE.begin(), reverseE.end());
  E += reverseE;

  assert(!ominoChecker.has_half_turn_tiling(B + C + D + E).empty());

  assert(ominoChecker.has_half_turn_tiling("NENNESSSEESWWWNW").empty());

  assert(ominoChecker.has_half_turn_tiling("WWWWWNNESEEENESS").empty());
}

void test__has_quarter_turn_tiling(){
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_quarter_turn_tiling(create_square(i)).empty());
  }

  assert(!ominoChecker.has_quarter_turn_tiling(create_long_rectangle(2)).empty());
  for (int i = 3; i <= 5; ++i) {
    assert(ominoChecker.has_quarter_turn_tiling(create_long_rectangle(i)).empty());
  }

  assert(!ominoChecker.has_quarter_turn_tiling("NNNESESSWW").empty());
  assert(ominoChecker.has_quarter_turn_tiling("NNENESESSWNWSW").empty());
}

void test__has_type_1_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_1_reflection_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_1_reflection_tiling(create_2_by_i_rectangle(i)).empty());
  }

  for (int i = 4; i <= 7; ++i) {
    assert(!ominoChecker.has_type_1_reflection_tiling(create_rectangle(3, i)).empty());
  }

  assert(ominoChecker.has_type_1_reflection_tiling("NWNENENESESESWSWNWSW").empty());
  assert(!ominoChecker.has_type_1_reflection_tiling("NNWNENNNENWNEEEENESEEESESESESWSWSWWWWNWSWWWW").empty());
}

void test__has_type_2_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_2_reflection_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_2_reflection_tiling(create_rectangle(2, i)).empty());
  }

  for (int i = 4; i <= 7; ++i) {
    assert(!ominoChecker.has_type_2_reflection_tiling(create_rectangle(3, i)).empty());
  }

  // tetris pieces (non-empty A, A hat)
  assert(!ominoChecker.has_type_2_reflection_tiling("NNENWNNNEENNEEEENEESESSWWSSSESWSSWNWWSWWWW").empty());
  assert(ominoChecker.has_type_2_reflection_tiling("NNENWNNNEENNEEEENEESESSWWSSSWSESSWNWWSWWWW").empty());

  // tetris pieces (empty A, A hat)
  assert(!ominoChecker.has_type_2_reflection_tiling("WNWNNNEESSEEWSWS").empty());
  assert(ominoChecker.has_type_2_reflection_tiling("WNWNNNEESSEEWWSS").empty());
}

void test__has_type_1_half_turn_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_1_half_turn_reflection_tiling(create_square(i)).empty());
  }

  // tetris pieces (non-empty A, A hat)
  assert(ominoChecker.has_type_1_half_turn_reflection_tiling("NNWNENNNENWNEEEENESEEESSSSSESWSSWWWNWSWWWW").empty());
  assert(!ominoChecker.has_type_1_half_turn_reflection_tiling("NNNWNENNNNENWNEEEENESEEESSSSSESWWSESWWWNWSWWWW").empty());

  // tetris pieces (empty A, A hat)
  assert(!ominoChecker.has_type_1_half_turn_reflection_tiling("NNNNNNESESESWSWWSW").empty());
  assert(!ominoChecker.has_type_1_half_turn_reflection_tiling("NNWNEENWNNESESESWSWWSW").empty());
  assert(ominoChecker.has_type_1_half_turn_reflection_tiling("NNWNENNNESESESWSWWSW").empty());
  assert(ominoChecker.has_type_1_half_turn_reflection_tiling("NNWNEENWNNESSEESWSWWSW").empty());
}

void test__has_type_2_half_turn_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!ominoChecker.has_type_2_half_turn_reflection_tiling(create_square(i)).empty());
  }
  // Tetris pieces (non-empty D, non-empty B, refl(B), empty A, C)
  assert(!ominoChecker.has_type_2_half_turn_reflection_tiling("NENWNENESESESWSWNWSW").empty());
  assert(ominoChecker.has_type_2_half_turn_reflection_tiling("NENWNENESESWSESWNWSW").empty());

  // Tetris pieces (non-empty D, non-empty B, refl(B), non-empty A, C)
  assert(!ominoChecker.has_type_2_half_turn_reflection_tiling("NENWWNENNENWNENESESESWSSSSSWNWSW").empty());
  assert(ominoChecker.has_type_2_half_turn_reflection_tiling("NENWWNENNENWNENESESESWSSESWSSWNWSW").empty());

  // Tetris pieces (empty D, non-empty B, refl(B) non-empty A, C)
  assert(!ominoChecker.has_type_2_half_turn_reflection_tiling("NENWWNENENESESSSSWNWSW").empty());
  assert(ominoChecker.has_type_2_half_turn_reflection_tiling("NENWWNENENESESSESWSWNWSW").empty());
}

int main() {
  // Test longest longest_match
  string S1 = "abcdef";
  string S2 = "abcfged";

  test__admissible_rotadrome_factors();
  test__admissible_reflect_square_factors();
  test__admissible_gapped_reflect_square_factor_pairs();

  // Test admissible_mirror_factors
  /*string boundary = "EESSENESSESWWWSWSWNNWWNEEENWNN";
  for (auto& p: admissible_mirror_factors(boundary)) {
    cout << p.first << " " << p.second << "\n";
  }*/

  test__has_translation_tiling();
  test__has_half_turn_tiling();
  test__has_quarter_turn_tiling();
  test__has_type_1_reflection_tiling();
  test__has_type_2_reflection_tiling();
  test__has_type_1_half_turn_reflection_tiling();
  test__has_type_2_half_turn_reflection_tiling();

 }
