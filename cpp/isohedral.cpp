#include "isohedral.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
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

template <typename T>
void test_equal(const string& test_name, T& actual, T& expected) {
    if (actual == expected) {
        cout << "[PASS] " << test_name << ": Actual = " << actual << ", Expected = " << expected << "\n";
    } else {
        cout << "[FAIL] " << test_name << ": Actual = " << actual << ", Expected = " << expected << "\n";
    }
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

void test__admissible_rotadrome_factors() {
  vector<Factor> factors = admissible_rotadrome_factors("NESW", 180);
  set<Factor> factorSet(factors.begin(), factors.end());
  set<Factor> wantSet = {{0, 0}, {1, 1}, {2, 2}, {3, 3}};
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
  /*for (int i = 1; i <= 5; ++i) {
    assert(!has_half_turn_tiling(create_square(i)).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    assert(!has_half_turn_tiling(create_long_rectangle(i)).empty());
    assert(!has_half_turn_tiling(create_tall_rectangle(i)).empty());
  }

  assert(!has_half_turn_tiling("NNESESSWWW").empty());
  assert(!has_half_turn_tiling("NNNESESSWW").empty());*/

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

int main() {
  // Test longest longest_match
  string S1 = "abcdef";
  string S2 = "abcfged";

  assert(longest_match(S1, S2, 5) == 3);
  assert(longest_match(S1, S2, 2) == 2);

  test__admissible_rotadrome_factors();

  // Test admissible_mirror_factors
  /*string boundary = "EESSENESSESWWWSWSWNNWWNEEENWNN";
  for (auto& p: admissible_mirror_factors(boundary)) {
    cout << p.first << " " << p.second << "\n";
  }*/

  test__has_translation_tiling();
  test__has_half_turn_tiling();

 }
