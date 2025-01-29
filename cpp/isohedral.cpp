#include "isohedral.h"

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

using namespace std;

using Factor = pair<int, int>;

void printFactor(const Factor& f) {
  cout << f.first << " " << f.second << "\n";
}

unordered_map<char, char> COMPLEMENT = {
  {'N', 'S'},
  {'S', 'N'},
  {'E', 'W'},
  {'W', 'E'}
};

unordered_map<char, char> CCW = {
  {'N', 'W'},
  {'E', 'N'},
  {'S', 'E'},
  {'W', 'S'}
};

unordered_map<int, unordered_map<char, char>> REFL = {
  {-45, {{'N', 'W'}, {'E', 'S'}, {'S', 'E'}, {'W', 'N'}}},
  {0, {{'N', 'S'}, {'E', 'E'}, {'S', 'N'}, {'W', 'W'}}},
  {45, {{'N', 'E'}, {'E', 'N'}, {'S', 'W'}, {'W', 'S'}}},
  {90, {{'N', 'N'}, {'E', 'W'}, {'S', 'S'}, {'W', 'E'}}},
};

int longest_match(const string& S1, const string& S2, size_t ub) {
  int i = 0;
  while (i < min({S1.length(), S2.length(), ub}) && S1[i] == S2[i]) {
    ++i;
  }
  return i;
}

vector<Factor> is_double_palindrome(const Factor& F, const vector<vector<Factor>>& palindrome_factor_starts, const vector<vector<Factor>>& palindrome_factor_ends, int n) {
  int F_len = F.second - F.first + 1 + n * (F.second < F.first);
  for (const auto& F1: palindrome_factor_starts[F.first]) {
    int F1_len = F1.second - F1.first + 1 + n * (F1.second < F1.first);
    if (F1_len == F_len) {
      return {F1};
    }
    for (const auto& F2: palindrome_factor_ends[F.second]) {
      int F2_len = F2.second - F2.first + 1 + n * (F2.second < F2.first);
      if (F_len == F1_len + F2_len) {
        return {F1, F2};
      }
    }
  }
  return {};
}

bool is_reflect_square_factor(const string& P, int i, int j, int theta) {
  int n = P.length();
  int l = j - i + 1 + n * (j < i);
  if (l % 2 != 0) {
    return false;
  }
  l /= 2;
  string reflected;
  /*cout << "i is: " << i << " l is: " << l << " sum is: " << i + l << "\n";
  cout << "Substring is: " << P.substr(i + l) << "\n";*/

  for (auto& c: P.substr(min(i + l, n))) {
    reflected += REFL[theta][c];
  }
  for (auto& c: P) {
    reflected += REFL[theta][c];
  }
  return (l == longest_match(P.substr(i) + P, reflected, l+1));
}

string inv_comp(const std::string& S) {
  std::string result;
  result.reserve(S.length());

  for (auto it = S.rbegin(); it != S.rend(); ++it) {
    result.push_back(COMPLEMENT[*it]);
  }
  return result;
}

vector<Factor> admissible_mirror_factors(const std::string& P) {
  int n = P.length();
  vector<Factor> factors;

  // Compute admissible mirror factors starting between letter pairs
  for (int i = 0; i < n; ++i) {
    int l = longest_match(inv_comp(P + P.substr(0, i)), P.substr((i + n / 2) % n) + P, n/4);
    int r = longest_match(P.substr(i) + P, inv_comp(P + P.substr(0, (i + n / 2) % n)), n/4);
    if (l == r && r > 0) {
      auto start = (i - l + n) % n;
      auto end = (i - 1 + r + n) % n;
      factors.push_back(make_pair(start, end));
    }
  }

  // Compute admissible mirror factors starting in middle of a letter
  for (int i = 0; i < n; ++i) {
    if (P[i] == COMPLEMENT[P[(i + n/2) % n]]) {
      int l = longest_match(inv_comp(P + P.substr(0, i)), P.substr((i + n/2 + 1) % n) + P, (n - 2)/4);
      int r = longest_match(P.substr((i+1)%n) + P, inv_comp(P + P.substr(0, (i + n/2) %n)), (n - 2)/4);
      if (l == r) {
        int start = (i - l + n) %n;
        int end = (i + r) % n;
        factors.push_back(make_pair(start, end));
      }
    }
  }
  return factors;
}

vector<pair<Factor, Factor>> admissible_gapped_mirror_factor_pairs(const string& P) {
  int n = P.length();
  vector<pair<Factor, Factor>> factor_pairs;
  // Compute admissible mirror factors starting between letter pairs
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int l = longest_match(inv_comp(P + P.substr(0, i)), P.substr(j) + P, (i + n - j) / 2);
      int r = longest_match(P.substr(i) + P, inv_comp(P + P.substr(0, j)), (j - i) / 2);

      if (l == r && r > 0) {
        factor_pairs.push_back({make_pair((i-l+n)%n, (i-1+r+n)%n), make_pair((j-l+n)%n, (j-1+r+n)%n)});
      }
    }
  }

  // Compute admissible mirror factors starting in the middle of a letter
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (P[i] == COMPLEMENT[P[j]]) {
        int l = longest_match(inv_comp(P + P.substr(0, i)), P.substr((j + 1) % n), (i + n - j - 1)/2);
        int r = longest_match(P.substr((i+1)%n) + P, inv_comp(P + P.substr(0, j)), (j-i-1)/2);
        if (l == r) {
          factor_pairs.push_back({make_pair((i-l+n)%n, (i+r)%n), make_pair((j-l+n)%n, (j+r)%n)});
        }
      }
    }
  }
  return factor_pairs;
}

vector<Factor> admissible_rotadrome_factors(const string& P, int theta) {
  int n = P.length();
  vector<Factor> factors;
  // Compute admissible palindrome factors starting between letter pairs
  for (int i = 0; i < n; ++i) {
    string firstString = P + P.substr(0, i);
    reverse(firstString.begin(), firstString.end());
    string rotString(n + (n - i), 'X');
    string secondString = P.substr(i) + P;
    for (int i = 0; i < secondString.length(); ++i) {
      rotString[i] = (theta == 90) ? CCW[secondString[i]] : secondString[i];
    }
    int l = longest_match(firstString, rotString, n/2);
    if (l > 0) {
      factors.push_back(make_pair((i-l+n)%n, (i-1+l)%n));
    }
  }
  // Compute admissible palindrome factors starting in middle of a letter
  if (theta == 180) {
    for (int i = 0; i < n; ++i) {
      string firstString = P + P.substr(0, i);
      reverse(firstString.begin(), firstString.end());
      int l = longest_match(firstString, P.substr((i+1)%n) + P, n/2);
      factors.push_back(make_pair((i-l+n)%n, (i+l)%n));
    }
  }
  return factors;
}

vector<Factor> admissible_reflect_square_factors(const string& P) {
  int n = P.length();

  set<Factor> factors;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      for (int theta = -45; theta <= 90; theta += 45) {
        if (is_reflect_square_factor(P, i, j, theta)) {
          factors.insert({i, j});
        }
      }
    }
  }
  return vector<Factor>(factors.begin(), factors.end());
}

vector<pair<Factor, Factor>> admissible_gapped_reflect_square_factor_pairs(const string& P, int theta) {
  int n = P.length();
  vector<pair<Factor, Factor>> factor_pairs;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      int d = min(j - i + n * (j < i), i - j + n * (i < j));
      string reflected;
      for (auto& c: P.substr(j)) {
        reflected += REFL[theta][c];
      }
      for (auto& c: P) {
        reflected += REFL[theta][c];
      }
      int l = longest_match(P.substr(i) + P, reflected, d+1);
      if (1 <= l && l <= d) {
        factor_pairs.push_back({{i, (i+l-1+n)%n}, {j, (j+l-1+n)%n}});
      }
    }
  }
  return factor_pairs;
}

vector<Factor> has_translation_tiling(const string& P) {
  int n = P.length();
  vector<Factor> factors = admissible_mirror_factors(P);
  vector<set<Factor>> factor_starts(n);
  vector<set<Factor>> factor_ends(n);
  for (auto& f: factors) {
    factor_starts[f.first].insert(f);
    factor_ends[f.second].insert(f);
  }

  for (int i = 0; i < n; ++i) {
    for (auto& A: factor_starts[i]) {
      for (auto& B: factor_starts[(A.second + 1) %n]) {
        int AB_len = A.second - A.first + 1 + n * (A.second < A.first) + B.second - B.first + 1 + n * (B.second < B.first);
        if (AB_len > n / 2) {
          continue;
        }
        if (AB_len == n / 2) {
          return {A, B};
        }
        Factor C = {(B.second + 1)%n, (A.first + n/2 - 1 + n)%n};
        if (factor_starts[C.first].find(C) != factor_starts[C.first].end()) {
          return {A, B, C};
        }
      }
    }
  }
  return {};
}

vector<Factor> has_half_turn_tiling(const string& P) {
  int n = P.length();
  auto mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P);
  /*cout << "Printing mirrors";
  for (auto& p: mirror_factor_pairs) {
    cout << "(";
    printFactor(p.first);
    cout << ")";
    cout << " (";
    printFactor(p.second);
    cout << ")\n";
  }*/
  vector<vector<Factor>> palindrome_factor_starts(n);
  vector<vector<Factor>> palindrome_factor_ends(n);
  for (auto& f: admissible_rotadrome_factors(P, 180)) {
    palindrome_factor_starts[f.first].push_back(f);
    palindrome_factor_ends[f.second].push_back(f);
  }

  Factor last_A;
  Factor last_A_hat;
  for (auto& p: mirror_factor_pairs)  {
    Factor A = p.first;
    Factor A_hat = p.second;

    Factor dpi1 = make_pair((A.second+1)%n, (A_hat.first-1+n)%n);
    Factor dpi2 = make_pair((A_hat.second+1)%n, (A.first-1+n)%n);
    vector<Factor> dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
    vector<Factor> dp2 = is_double_palindrome(dpi2, palindrome_factor_starts, palindrome_factor_ends, n);

    if (!dp1.empty() && !dp2.empty()) {
      vector<Factor> result;
      result.push_back(A);
      result.insert(result.end(), dp1.begin(), dp1.end());
      result.push_back(A_hat);
      result.insert(result.end(), dp2.begin(), dp2.end());
      return result;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      Factor dpi1 = make_pair(i, (j-1+n)%n);
      Factor dpi2 = make_pair(j%n, (i-1+n)%n);
      vector<Factor> dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
      vector<Factor> dp2 = is_double_palindrome(dpi2, palindrome_factor_starts, palindrome_factor_ends, n);
      if (!dp1.empty() && !dp2.empty()) {
        vector<Factor> result;
        result.push_back(last_A);
        result.insert(result.end(), dp1.begin(), dp1.end());
        result.push_back(last_A_hat);
        result.insert(result.end(), dp2.begin(), dp2.end());
        return result;
      }
    }
  }
  return {};
}

vector<Factor> has_quarter_turn_tiling(const string& P) {
  int n = P.length();
  vector<Factor> palin_factors = admissible_rotadrome_factors(P, 180);
  vector<Factor> ninety_factors = admissible_rotadrome_factors(P, 90);
  vector<vector<Factor>> ninety_factor_starts(n);
  for (auto& f: ninety_factors) {
    ninety_factor_starts[f.first].push_back(f);
  }
  // Factorizations with non-empty palindrome factor
  for (auto& C: palin_factors) {
    int C_len = C.second - C.first + 1 + n * (C.second < C.first);
    for (auto& A: ninety_factor_starts[(C.second + 1)%n]) {
      int A_len = A.second - A.first + 1 + n * (A.second < A.first);
      // One empty 90-drome factor
      if (A_len + C_len == n) {
        return {C, A};
      }
      // No empty 90-drome factors
      for (auto& B: ninety_factor_starts[(A.second+1)%n]) {
        int B_len = B.second - B.first + 1 + n * (B.second < B.first);
        if (A_len + B_len + C_len == n) {
          return {C, A, B};
        }
      }
    }
  }
  // Factorizations with empty palindrome factor
  for (auto& A: ninety_factors) {
    for (auto& B: ninety_factor_starts[(A.second + 1)%n]) {
      int AB_len = A.second - A.first + 1 + n * (A.second < A.first) + B.second - B.first + 1 + n * (B.second < B.first);
      if (AB_len == n) {
        return {A, B};
      }
    }
  }
  return {};
}

vector<Factor> has_type_1_reflection_tiling(const string& P) {
  int n = P.length();
  vector<pair<Factor, Factor>> mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P);
  vector<Factor> reflect_square_factors = admissible_reflect_square_factors(P);
  vector<set<Factor>> reflect_square_factor_starts(n);
  for (auto& f: reflect_square_factors){
    reflect_square_factor_starts[f.first].insert(f);
  }
  for (auto& p: mirror_factor_pairs) {
    Factor A = p.first;
    Factor A_hat = p.second;
    Factor rem1 = {(A.second+1)%n, (A_hat.first-1+n)%n};
    Factor rem2 = {(A_hat.second+1)%n, (A.first-1+n)%n};
    if (reflect_square_factor_starts[rem1.first].find(rem1) != reflect_square_factor_starts[rem1.first].end() && reflect_square_factor_starts[rem2.first].find(rem2) != reflect_square_factor_starts[rem2.first].end()) {
      return {A, rem1, A_hat, rem2};
    }
  }
  for (auto& f: reflect_square_factors) {
    int f_len = f.second - f.first + 1 + n * (f.second < f.first);
    for (auto& of: reflect_square_factor_starts[(f.second+1)%n]) {
      int of_len = of.second - of.first + 1 + n * (of.second < of.first);
      if (f_len + of_len == n) {
        return {f, of};
      }
    }
  }
  return {};
}

vector<Factor> has_type_2_reflection_tiling(const string& P) {
  int n = P.length();
  vector<Factor> mirror_factors = admissible_mirror_factors(P);
  for (int theta = -45; theta <= 90; theta += 45) {
    map<Factor, vector<pair<Factor, Factor>>> reflect_factor_tips;
    for (auto& p: admissible_gapped_reflect_square_factor_pairs(P, theta)) {
      Factor f = p.first;
      Factor cf = p.second;
      reflect_factor_tips[make_pair(f.first, cf.second)].push_back({f, cf});
      reflect_factor_tips[make_pair(cf.first, f.second)].push_back({cf, f});
    }
    for (auto& A: mirror_factors) {
      int A_len = A.second - A.first + 1 + n * (A.second < A.first);
      Factor rem1 = {(A.second+1)%n, (A.first-1+n)%n};
      for (auto& p: reflect_factor_tips[rem1]) {
        Factor f1 = p.first;
        Factor cf1 = p.second;
        Factor rem2 = {(f1.second+1)%n, (cf1.first-1+n)%n};
        int rem2_len = rem2.second - rem2.first + 1 + n * (rem2.second < rem2.first);
        if (rem2_len == A_len) {
          return {A, f1, rem2, cf1};
        }
        for (auto& p2: reflect_factor_tips[rem2]) {
          Factor f2 = p2.first;
          Factor cf2 = p2.second;
          Factor rem3 = {(f2.second+1)%n, (cf2.first-1+n)%n};
          int rem3_len = rem3.second - rem3.first + 1 + n * (rem3.second < rem3.first);
          if (rem3_len == A_len) {
            return {A, f1, f2, rem3, cf2, cf1};
          }
        }
      }
    }
    for (int i = 0; i < n; ++i) {
      for (auto& p: reflect_factor_tips[{(i+1)%n, i}]) {
        Factor f1 = p.first;
        Factor cf1 = p.second;
        Factor rem2 = {(f1.second+1)%n, (cf1.first-1+n)%n};
        for (auto& p2: reflect_factor_tips[rem2]) {
          Factor f2 = p2.first;
          Factor cf2 = p2.second;
          if ((f2.second+1)%n == cf2.first) {
            return {f1, f2, cf2, cf1};
          }
        }
      }
    }
  }
  return {};
}

vector<Factor> has_type_1_half_turn_reflection_tiling(const string& P) {
  int n = P.length();
  vector<pair<Factor, Factor>> mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P);
  vector<Factor> palin_factors = admissible_rotadrome_factors(P, 180);
  vector<Factor> reflect_square_factors = admissible_reflect_square_factors(P);

  vector<vector<Factor>> palindrome_factor_starts(n);
  vector<vector<Factor>> palindrome_factor_ends(n);
  for (auto& f: palin_factors) {
    palindrome_factor_starts[f.first].push_back(f);
    palindrome_factor_ends[f.second].push_back(f);
  }
  vector<vector<Factor>> reflect_square_starts(n);
  for (auto& f: reflect_square_factors) {
    reflect_square_starts[f.first].push_back(f);
  }

  // Non-empty mirrors 
  for (auto& p: mirror_factor_pairs) {
    Factor A = p.first;
    Factor A_hat = p.second;
    Factor dpi1 = {(A.second+1)%n, (A_hat.first-1+n)%n};
    Factor dpi2 = {(A_hat.second+1)%n, (A.first-1+n)%n};
    vector<Factor> dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
    // Other variant (dpi2 is the double palindrome) checked later by symmetry of A, A_hat
    vector<Factor>& starts = reflect_square_starts[dpi2.first];
    if (!dp1.empty() && (find(starts.begin(), starts.end(), dpi2) != starts.end())) {
      vector<Factor> result;
      result.push_back(A);
      result.insert(result.end(), dp1.begin(), dp1.end());
      result.push_back(A_hat);
      result.push_back(dpi2);
      return result;
    }
  }
  
  // Empty mirrors
  for (auto& dp1: reflect_square_factors) {
    vector<Factor> dp2 = is_double_palindrome({(dp1.second+1)%n, (dp1.first-1+n)%n}, palindrome_factor_starts, palindrome_factor_ends, n);
    if (!dp2.empty()) {
      vector<Factor> result = {dp1};
      result.insert(result.end(), dp2.begin(), dp2.end());
      return result;
    }
  }
  return {};
}

/****************************************************
* Code for tests 
* This section contains:
*   - Helper functions 
*   - Functions that test functionality.
*   - A main function which runs our test suit. 
****************************************************/

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
  vector<Factor> factors = admissible_rotadrome_factors("NESW", 180);
  set<Factor> factorSet(factors.begin(), factors.end());
  set<Factor> wantSet = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
  assert(factorSet == wantSet);
}

void test__admissible_reflect_square_factors() {
  vector<Factor> factors = admissible_reflect_square_factors("NEESWW");
  set<Factor> factorSet(factors.begin(), factors.end());
  set<Factor> wantSet = {{0, 1}, {2, 3}, {3, 4}, {5, 0}, {1, 2}, {4, 5}};
  assert(factorSet == wantSet);

  factors = admissible_reflect_square_factors("NNEESSWW");
  factorSet = set<Factor>(factors.begin(), factors.end());
  wantSet = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 0}, {0, 3}, {2, 5}, {4, 7}, {6, 1}};
  assert(factorSet == wantSet);
}

void test__has_translation_tiling() {
  // Test has_translation_tiling
  for (int i = 1; i <= 5; ++i) {
    assert(!has_translation_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!has_translation_tiling(create_long_rectangle(i)).empty());

    // Rectangle in other direction 
    assert(!has_translation_tiling(create_tall_rectangle(i)).empty());
  }

  assert(!has_translation_tiling("NNESESSWNW").empty());

  assert(!has_translation_tiling("NNNESESSWW").empty());
  assert(!has_translation_tiling("NENENESEESWSWNWSWW").empty());
  assert(has_translation_tiling("NENENESEESSWWNWSWW").empty());

  assert(has_translation_tiling("NENNESSSEESWWWNW").empty());
  assert(has_translation_tiling("WWWNNESENESS").empty());
}

void test__has_half_turn_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!has_half_turn_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!has_half_turn_tiling(create_long_rectangle(i)).empty());
    assert(!has_half_turn_tiling(create_tall_rectangle(i)).empty());
  }

  assert(!has_half_turn_tiling("NNESESSWWW").empty());
  assert(!has_half_turn_tiling("NNNESESSWW").empty());

  assert(!has_half_turn_tiling("NWNEENWNENESESESSWSWNWSW").empty());

  assert(!has_half_turn_tiling("NENENESEESWSWNWSWW").empty());
  assert(!has_half_turn_tiling("WWWNNESENESS").empty());

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

  assert(!has_half_turn_tiling(B + C + D + E).empty());

  assert(has_half_turn_tiling("NENNESSSEESWWWNW").empty());

  assert(has_half_turn_tiling("WWWWWNNESEEENESS").empty());
}

void test__has_quarter_turn_tiling(){
  for (int i = 1; i <= 5; ++i) {
    assert(!has_quarter_turn_tiling(create_square(i)).empty());
  }

  assert(!has_quarter_turn_tiling(create_long_rectangle(2)).empty());
  for (int i = 3; i <= 5; ++i) {
    assert(has_quarter_turn_tiling(create_long_rectangle(i)).empty());
  }

  assert(!has_quarter_turn_tiling("NNNESESSWW").empty());
  assert(has_quarter_turn_tiling("NNENESESSWNWSW").empty());
}

void test__has_type_1_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!has_type_1_reflection_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!has_type_1_reflection_tiling(create_2_by_i_rectangle(i)).empty());
  }

  for (int i = 4; i <= 7; ++i) {
    assert(!has_type_1_reflection_tiling(create_rectangle(3, i)).empty());
  }

  assert(has_type_1_reflection_tiling("NWNENENESESESWSWNWSW").empty());
  assert(!has_type_1_reflection_tiling("NNWNENNNENWNEEEENESEEESESESESWSWSWWWWNWSWWWW").empty());
}

void test__has_type_2_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!has_type_2_reflection_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!has_type_2_reflection_tiling(create_rectangle(2, i)).empty());
  }

  for (int i = 4; i <= 7; ++i) {
    assert(!has_type_2_reflection_tiling(create_rectangle(3, i)).empty());
  }

  // tetris pieces (non-empty A, A hat)
  assert(!has_type_2_reflection_tiling("NNENWNNNEENNEEEENEESESSWWSSSESWSSWNWWSWWWW").empty());
  assert(has_type_2_reflection_tiling("NNENWNNNEENNEEEENEESESSWWSSSWSESSWNWWSWWWW").empty());

  // tetris pieces (empty A, A hat)
  assert(!has_type_2_reflection_tiling("WNWNNNEESSEEWSWS").empty());
  assert(has_type_2_reflection_tiling("WNWNNNEESSEEWWSS").empty());
}

void test__has_type_1_half_turn_reflection_tiling() {
  for (int i = 1; i <= 5; ++i) {
    assert(!has_type_1_half_turn_reflection_tiling(create_square(i)).empty());
  }

  // tetris pieces (non-empty A, A hat)
  assert(has_type_1_half_turn_reflection_tiling("NNWNENNNENWNEEEENESEEESSSSSESWSSWWWNWSWWWW").empty());
  assert(!has_type_1_half_turn_reflection_tiling("NNNWNENNNNENWNEEEENESEEESSSSSESWWSESWWWNWSWWWW").empty());

  // tetris pieces (empty A, A hat)
  assert(!has_type_1_half_turn_reflection_tiling("NNNNNNESESESWSWWSW").empty());
  assert(!has_type_1_half_turn_reflection_tiling("NNWNEENWNNESESESWSWWSW").empty());
  assert(has_type_1_half_turn_reflection_tiling("NNWNENNNESESESWSWWSW").empty());
  assert(has_type_1_half_turn_reflection_tiling("NNWNEENWNNESSEESWSWWSW").empty());

}

int main() {
  // Test longest longest_match
  string S1 = "abcdef";
  string S2 = "abcfged";

  assert(longest_match(S1, S2, 5) == 3);
  assert(longest_match(S1, S2, 2) == 2);

  test__admissible_rotadrome_factors();
  test__admissible_reflect_square_factors();

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

 }
