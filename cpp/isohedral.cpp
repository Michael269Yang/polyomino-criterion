#include "isohedral.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

unordered_map<char, char> COMPLEMENT = {
  {'N', 'S'},
  {'S', 'N'},
  {'E', 'W'},
  {'W', 'E'}
};

int longest_match(const string& S1, const string& S2, size_t ub) {
  int i = 0;
  while (i < min({S1.length(), S2.length(), ub}) && S1[i] == S2[i]) {
    ++i;
  }
  return i;
}

string inv_comp(const std::string& S) {
  std::string result;
  result.reserve(S.length());

  for (auto it = S.rbegin(); it != S.rend(); ++it) {
    result.push_back(COMPLEMENT[*it]);
  }
  return result;
}

vector<pair<int, int>> admissible_mirror_factors(const std::string& P) {
  int n = P.length();
  vector<pair<int, int>> factors;

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

vector<pair<int, int>> has_translation_tiling(const string& P) {
  int n = P.length();
  vector<pair<int, int>> factors = admissible_mirror_factors(P);
  vector<set<pair<int, int>>> factor_starts(n);
  vector<set<pair<int, int>>> factor_ends(n);
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
        pair<int, int> C = {(B.second + 1)%n, (A.first + n/2 - 1 + n)%n};
        if (factor_starts[C.first].find(C) != factor_starts[C.first].end()) {
          return {A, B, C};
        }
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


int main() {
  // Test longest longest_match
  string S1 = "abcdef";
  string S2 = "abcfged";

  assert(longest_match(S1, S2, 5) == 3);
  assert(longest_match(S1, S2, 2) == 2);

  // Test admissible_mirror_factors
  string boundary = "EESSENESSESWWWSWSWNNWWNEEENWNN";
  for (auto& p: admissible_mirror_factors(boundary)) {
    cout << p.first << " " << p.second << "\n";
  }

  // Test has_translation_tiling
  for (int i = 1; i <= 5; ++i) {
    string square(4 * i, ' ');
    fill(square.begin(), square.begin() + i, 'N');
    fill(square.begin() + i, square.begin() + 2 * i, 'E');
    fill(square.begin() + 2 * i, square.begin() + 3 * i, 'S');
    fill(square.begin() + 3 * i, square.end(), 'W');
    assert(!has_translation_tiling(square).empty());
  }

  for (int i = 1; i <= 5; ++i) {
    string rectangle(2 * i + 2, ' ');
    rectangle[0] = 'N';
    fill(rectangle.begin() + 1, rectangle.begin() + i + 1, 'E');
    rectangle[i + 1] = 'S';
    fill(rectangle.begin() + i + 2, rectangle.end(), 'W');
    assert(!has_translation_tiling(rectangle).empty());

    // Rectangle in other direction 
    fill(rectangle.begin(), rectangle.begin() + i, 'N');
    rectangle[i] = 'E';
    fill(rectangle.begin() + i + 1, rectangle.begin() + 2 * i + 1, 'S');
    rectangle[2 * i + 1] = 'W';
    assert(!has_translation_tiling(rectangle).empty());
  }

  assert(!has_translation_tiling("NNESESSWNW").empty());

  assert(!has_translation_tiling("NNNESESSWW").empty());
  assert(!has_translation_tiling("NENENESEESWSWNWSWW").empty());
  assert(has_translation_tiling("NENENESEESSWWNWSWW").empty());

  assert(has_translation_tiling("NENNESSSEESWWWNW").empty());
  assert(has_translation_tiling("WWWNNESENESS").empty());
}
