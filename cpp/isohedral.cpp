#include "boundary.h"
#include "isohedral.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

void printFactor(const Factor& f) {
  cout << f.first << " " << f.second << "\n";
}

int longest_match(const boundaryword& S1, const boundaryword& S2, size_t ub) {
  int i = 0;
  while (i < min({S1.size(), S2.size(), ub}) && S1[i] == S2[i]) {
    ++i;
  }
  return i;
}

bool is_double_palindrome(const Factor& F, const vector<vector<Factor>>& palindrome_factor_starts, const vector<vector<Factor>>& palindrome_factor_ends, int n) {
  int F_len = F.second - F.first + 1 + n * (F.second < F.first);
  for (const auto& F1: palindrome_factor_starts[F.first]) {
    int F1_len = F1.second - F1.first + 1 + n * (F1.second < F1.first);
    if (F1_len == F_len) {
      return true;
    }
    for (const auto& F2: palindrome_factor_ends[F.second]) {
      int F2_len = F2.second - F2.first + 1 + n * (F2.second < F2.first);
      if (F_len == F1_len + F2_len) {
        return true;
      }
    }
  }
  return false;
}

bool IsohedralChecker::is_reflect_square_factor(const boundaryword& P, int i, int j, int theta) {
  int n = P.size();
  int l = j - i + 1 + n * (j < i);
  if (l % 2 != 0) {
    return false;
  }
  l /= 2;
  boundaryword reflected;

  for (auto& c: boundaryword(P.begin() + ((i+l)%n), P.end())) {
    reflected.push_back(REFL[theta][c]);
  }
  for (auto& c: P) {
    reflected.push_back(REFL[theta][c]);
  }
  return (l == longest_match(slice(P, i) + P, reflected, l+1));
}

boundaryword IsohedralChecker::inv_comp(const boundaryword& S) {
  boundaryword result;
  result.reserve(S.size());

  for (auto it = S.rbegin(); it != S.rend(); ++it) {
    result.push_back(COMPLEMENT[*it]);
  }
  return result;
}

pair<int, int> IsohedralChecker::iteratedCcw(pair<int, int> dir, int numIters) {
  if (numIters < 0) return {0, 0};
  for (int i = 0; i < numIters; ++i) {
    dir = CCW[dir];
  }
  return dir;
}

FactorArray<2*MAX_BND> IsohedralChecker::admissible_mirror_factors(const boundaryword& P) {
  int n = P.size();
  FactorArray<2*MAX_BND> factors;
  size_t factor_idx = 0;

  // Compute admissible mirror factors starting between letter pairs
  for (int i = 0; i < n; ++i) {
    int l = longest_match(inv_comp(P + slice(P, 0, i)), slice(P, (i + n / 2) %n) + P, n/4);
    int r = longest_match(slice(P, i) + P, inv_comp(P + slice(P, 0, (i + n/2) %n)), n/4);
    if (l == r && r > 0) {
      auto start = (i - l + n) % n;
      auto end = (i - 1 + r + n) % n;
      factors[factor_idx++] = make_pair(start, end);
    }
  }

  // Compute admissible mirror factors starting in middle of a letter
  for (int i = 0; i < n; ++i) {
    if (P[i] == COMPLEMENT[P[(i + n/2) % n]]) {
      int l = longest_match(inv_comp(P + slice(P, 0, i)), slice(P, (i + n/2 + 1) %n) + P, (n - 2)/4);
      int r = longest_match(slice(P, (i+1)%n) + P, inv_comp(P + slice(P, 0, (i + n/2)%n)), (n-2)/4);
      if (l == r) {
        int start = (i - l + n) %n;
        int end = (i + r) % n;
        factors[factor_idx++] = make_pair(start, end);
      }
    }
  }
  factors.filled_count = factor_idx;
  return factors;
}

FactorPairArray<MAX_BND*MAX_BND> IsohedralChecker::admissible_gapped_mirror_factor_pairs(const boundaryword& P) {
  int n = P.size();
  FactorPairArray<MAX_BND*MAX_BND> factor_pairs;
  size_t factor_idx = 0;
  // Compute admissible mirror factors starting between letter pairs
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int l = longest_match(inv_comp(P + slice(P, 0, i)), slice(P, j)+ P, (i + n - j) / 2);
      int r = longest_match(slice(P, i) + P, inv_comp(P + slice(P, 0, j)), (j - i) / 2);

      if (l == r && r > 0) {
        factor_pairs[factor_idx++] = {make_pair((i-l+n)%n, (i-1+r+n)%n), make_pair((j-l+n)%n, (j-1+r+n)%n)};
      }
    }
  }

  // Compute admissible mirror factors starting in the middle of a letter
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (P[i] == COMPLEMENT[P[j]]) {
        int l = longest_match(inv_comp(P + slice(P, 0, i)), slice(P, (j + 1) % n) + P, (i + n - j - 1)/2);
        int r = longest_match(slice(P, (i+1)%n) + P, inv_comp(P + slice(P, 0, j)), (j-i-1)/2);
        if (l == r) {
          factor_pairs[factor_idx++] = {make_pair((i-l+n)%n, (i+r)%n), make_pair((j-l+n)%n, (j+r)%n)};
        }
      }
    }
  }
  factor_pairs.filled_count = factor_idx;
  return factor_pairs;
}

/*string convertToStr(const boundaryword& P) {
  pair<int, int> E = {3, 0};
  pair<int, int> NE = {0, 3};
  pair<int, int> NW = {-3, 3};
  pair<int, int> W = {-3, 0};
  pair<int, int> SW = {0, -3};
  pair<int, int> SE = {3, -3};
  map<pair<int, int>, char> theMap = {
    {E, 'E'}, {NE, 'R'}, {NW, 'L'}, {W, 'W'}, {SW, 'l'}, {SE, 'r'}
  };

  string ret;
  for (auto& p: P) {
    ret += theMap[p];
  }
  return ret;
}*/


string convertToStr(const boundaryword& P) {
  pair<int, int> U = {-1, 2};
  pair<int, int> D = {1, -2};
  pair<int, int> R = {1, 1};
  pair<int, int> r = {2, -1};
  pair<int, int> L = {-2, 1};
  pair<int, int> l = {-1, -1};

  map<pair<int, int>, char> theMap = {
    {U, 'U'}, {D, 'D'}, {R, 'R'}, {r, 'r'}, {L, 'L'}, {l, 'l'}
  };

  string ret;
  for (auto& p: P) {
    ret += theMap[p];
  }
  return ret;
}

FactorArray<2*MAX_BND> IsohedralChecker::admissible_rotadrome_factors(const boundaryword& P, int theta) {
  int n = P.size();
  theta = 180 - theta;
  int numCcw = (theta % minAngle == 0) ? theta / minAngle : -1;
  FactorArray<2*MAX_BND> factors;
  size_t factor_idx = 0;

  for (int i = 0; i < n; ++i) {
    boundaryword firstString = P + slice(P, 0, i);
    reverse(firstString.begin(), firstString.end());
    boundaryword rotString(n + (n - i), {0, 0});
    boundaryword secondString = slice(P, i)+ P;
    for (int i = 0; i < secondString.size(); ++i) {
      rotString[i] = iteratedCcw(secondString[i], numCcw);
    }
    int l = longest_match(firstString, rotString, n/2);
    if (l > 0) {
      factors[factor_idx++] = make_pair((i-l+n)%n, (i-1+l)%n);
    }
  }
  // Compute admissible palindrome factors starting in middle of a letter
  if (theta == 0) {
    for (int i = 0; i < n; ++i) {
      boundaryword firstString = P + slice(P, 0, i);
      reverse(firstString.begin(), firstString.end());
      int l = longest_match(firstString, slice(P, (i+1)%n) + P, n/2);
      factors[factor_idx++] = make_pair((i-l+n)%n, (i+l)%n);
    }
  }
  factors.filled_count = factor_idx;
  return factors;
}

FactorArray<MAX_BND*MAX_BND> IsohedralChecker::admissible_reflect_square_factors(const boundaryword& P) {
  int n = P.size();

  FactorArray<MAX_BND*MAX_BND> factors;
  size_t factor_idx = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      for (auto &p: REFL) {
        int theta = p.first;
        if (is_reflect_square_factor(P, i, j, theta)) {
          factors[factor_idx++] = {i, j};
        }
      }
    }
  }
  factors.filled_count = factor_idx;
  return factors;
}

vector<pair<Factor, Factor>> IsohedralChecker::admissible_gapped_reflect_square_factor_pairs(const boundaryword& P, int theta) {
  int n = P.size();
  vector<pair<Factor, Factor>> factor_pairs;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int d = min(j - i + n * (j < i), i - j + n * (i < j));
      boundaryword reflected;
      for (auto& c: slice(P, j)) {
        reflected.push_back(REFL[theta][c]);
      }
      for (auto& c: P) {
        reflected.push_back(REFL[theta][c]);
      }
      int l = longest_match(slice(P, i) + P, reflected, d+1);
      if (1 <= l && l <= d) {
        factor_pairs.push_back({{i, (i+l-1+n)%n}, {j, (j+l-1+n)%n}});
      }
    }
  }
  return factor_pairs;
}

bool IsohedralChecker::has_translation_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& mirror_factors) {
  int n = P.size();
  vector<set<Factor>> factor_starts(n);
  vector<set<Factor>> factor_ends(n);
  for (const auto& f: mirror_factors) {
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
          return true;
        }
        Factor C = {(B.second + 1)%n, (A.first + n/2 - 1 + n)%n};
        if (factor_starts[C.first].find(C) != factor_starts[C.first].end()) {
          return true;
        }
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_half_turn_tiling(const boundaryword& P, const FactorPairArray<MAX_BND*MAX_BND>& mirror_factor_pairs, const FactorArray<2*MAX_BND>& palin_factors) {
  int n = P.size();
  vector<vector<Factor>> palindrome_factor_starts(n);
  vector<vector<Factor>> palindrome_factor_ends(n);
  for (auto& f: palin_factors) {
    palindrome_factor_starts[f.first].push_back(f);
    palindrome_factor_ends[f.second].push_back(f);
  }

  for (auto& p: mirror_factor_pairs)  {
    Factor A = p.first;
    Factor A_hat = p.second;

    Factor dpi1 = make_pair((A.second+1)%n, (A_hat.first-1+n)%n);
    Factor dpi2 = make_pair((A_hat.second+1)%n, (A.first-1+n)%n);
    bool dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
    bool dp2 = is_double_palindrome(dpi2, palindrome_factor_starts, palindrome_factor_ends, n);

    if (dp1 && dp2) {
      return true;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      Factor dpi1 = make_pair(i, (j-1+n)%n);
      Factor dpi2 = make_pair(j%n, (i-1+n)%n);
      bool dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
      bool dp2 = is_double_palindrome(dpi2, palindrome_factor_starts, palindrome_factor_ends, n);
      if (dp1 && dp2) {
        return true;
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_quarter_turn_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& ninety_factors, const FactorArray<2*MAX_BND>& palin_factors) {
  int n = P.size();
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
        return true;
      }
      // No empty 90-drome factors
      for (auto& B: ninety_factor_starts[(A.second+1)%n]) {
        int B_len = B.second - B.first + 1 + n * (B.second < B.first);
        if (A_len + B_len + C_len == n) {
          return true;
        }
      }
    }
  }
  // Factorizations with empty palindrome factor
  for (auto& A: ninety_factors) {
    for (auto& B: ninety_factor_starts[(A.second + 1)%n]) {
      int AB_len = A.second - A.first + 1 + n * (A.second < A.first) + B.second - B.first + 1 + n * (B.second < B.first);
      if (AB_len == n) {
        return true;
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_type_1_reflection_tiling(const boundaryword& P, const FactorArray<MAX_BND*MAX_BND>& reflect_square_factors, const FactorPairArray<MAX_BND*MAX_BND>& mirror_factor_pairs) {
  int n = P.size();
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
      return true;
    }
  }
  for (auto& f: reflect_square_factors) {
    int f_len = f.second - f.first + 1 + n * (f.second < f.first);
    for (auto& of: reflect_square_factor_starts[(f.second+1)%n]) {
      int of_len = of.second - of.first + 1 + n * (of.second < of.first);
      if (f_len + of_len == n) {
        return true;
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_type_2_reflection_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& mirror_factors) {
  int n = P.size();
  for (auto& p: REFL) {
    int theta = p.first;
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
          return true;
        }
        for (auto& p2: reflect_factor_tips[rem2]) {
          Factor f2 = p2.first;
          Factor cf2 = p2.second;
          Factor rem3 = {(f2.second+1)%n, (cf2.first-1+n)%n};
          int rem3_len = rem3.second - rem3.first + 1 + n * (rem3.second < rem3.first);
          if (rem3_len == A_len) {
            return true;
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
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_type_1_half_turn_reflection_tiling(const boundaryword& P, const FactorPairArray<MAX_BND*MAX_BND>& partial_mirror_factor_pairs, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<MAX_BND*MAX_BND>& reflect_square_factors) {
  int n = P.size();
  // Factorization A B C A_hat D f_theta(D) is not symmetric so we need both orderings of each pair.
  vector<pair<Factor, Factor>> mirror_factor_pairs;
  mirror_factor_pairs.reserve(2 * partial_mirror_factor_pairs.filled_count);
  for (auto& p: partial_mirror_factor_pairs) {
    mirror_factor_pairs.push_back(p);
    mirror_factor_pairs.push_back({p.second, p.first});
  }
  
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
    bool dp1 = is_double_palindrome(dpi1, palindrome_factor_starts, palindrome_factor_ends, n);
    // Other variant (dpi2 is the double palindrome) checked later by symmetry of A, A_hat
    vector<Factor>& starts = reflect_square_starts[dpi2.first];
    if (dp1 && (find(starts.begin(), starts.end(), dpi2) != starts.end())) {
      return true;
    }
  }
  
  // Empty mirrors
  for (auto& dp1: reflect_square_factors) {
    bool dp2 = is_double_palindrome({(dp1.second+1)%n, (dp1.first-1+n)%n}, palindrome_factor_starts, palindrome_factor_ends, n);
    if (dp2) {
      return true;
    }
  }
  return false;
}

bool IsohedralChecker::has_type_2_half_turn_reflection_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors) {
  int n = P.size();
  map<int, vector<pair<Factor, Factor>>> reflect_factor_pairs;
  for (auto& p: REFL) {
    int theta = p.first;
    reflect_factor_pairs[theta] = admissible_gapped_reflect_square_factor_pairs(P, theta);
    // Double up, both for tips and for f1, cf1 iterating since B, refl(B) are not symmetric
    vector<pair<Factor, Factor>> reversed_pairs;
    for (auto& fp: reflect_factor_pairs[theta]) {
      reversed_pairs.push_back({fp.second, fp.first});
    }
    reflect_factor_pairs[theta].insert(reflect_factor_pairs[theta].end(), reversed_pairs.begin(), reversed_pairs.end());
  }

  vector<pair<int, int>> theta_pairs;
  for (auto&p: REFL){
    int theta = p.first;
    int perpendicular = (theta > 0) ? theta - 90: theta + 90;
    theta_pairs.push_back({theta, perpendicular});
  }
  for (auto& p: theta_pairs) {
    int theta1 = p.first;
    int theta2 = p.second;
    map<Factor, vector<pair<Factor, Factor>>> reflect_factor_tips;
    for (auto& pf: reflect_factor_pairs[theta2]) {
      Factor f = pf.first;
      Factor cf = pf.second;
      reflect_factor_tips[{f.first, cf.second}].push_back(make_pair(f, cf));
    }
    // Non-empty B, refl(B)
    for (auto& pf: reflect_factor_pairs[theta1]) {
      Factor f1 = pf.first;
      Factor cf1 = pf.second;
      Factor dpi1 = {(f1.second+1)%n, (cf1.first-1+n)%n};
      Factor dpi2 = {(cf1.second+1)%n, (f1.first-1+n)%n};
      // Empty D, refl(D)
      if ((find(palin_factors.begin(), palin_factors.end(), dpi1) != palin_factors.end()) && (find(palin_factors.begin(), palin_factors.end(), dpi2) != palin_factors.end())) {
        return true;
      }
      // Non-empty D, refl(D)
      for (auto & pr: reflect_factor_tips[{dpi1.first, dpi2.second}]) {
        Factor f2 = pr.first;
        Factor cf2 = pr.second;
        vector<Factor> rem1f;
        if ((f2.second+1)%n == cf1.first) { // |A| = 0
          rem1f = {f2};
        }
        Factor rem1i = {(f2.second+1)%n, dpi1.second};

        if (rem1f.empty() && (find(palin_factors.begin(), palin_factors.end(), rem1i) != palin_factors.end())) {
          rem1f =  {f2, rem1i};
        }

        vector<Factor> rem2f;
        if (cf1.second == (cf2.first-1+n)%n) { // |C| = 0
          rem2f = {cf2};
        }
        Factor rem2i = {dpi2.first, (cf2.first-1+n)%n};
        if (rem2f.empty() && (find(palin_factors.begin(), palin_factors.end(), rem2i) != palin_factors.end())) {
          rem2f = {rem2i, cf2};
        }

        if (!rem1f.empty() && !rem2f.empty()) {
          return true;
        }
      }
      // Empty B, refl(B): D refl(D) A C with A, C, palindromes
      for (int i = 0; i < n; ++i) {
        for (auto& pr: reflect_factor_tips[{(i+1)%n, i}]){
          Factor f2 = pr.first;
          Factor cf2 = pr.second;
          if ((f2.second+1)%n == cf2.first && (cf2.second+1)%n == f2.first) {
            return true;
          }
          for (auto& p1: palin_factors) {
            if (p1.first != (cf2.second+1)%n) {
              continue;;
            }
            if ((p1.second+1)%n == f2.first) {
              return true;
            }
            for (auto& p2: palin_factors) {
              if (p2.first != (p1.second+1)%n || (p2.second+1)%n != f2.first) {
                continue;
              }
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_case_7_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& onetwenty_factors) {
  // A t_120(A) B t_120(B) C t_120(C)
  int n = P.size();
  vector<set<Factor>> factor_starts(n);
  for (auto& f: onetwenty_factors) {
    factor_starts[f.first].insert(f);
  }

  for (auto& A: onetwenty_factors) {
    int A_len = A.second - A.first + 1 + n * (A.second < A.first);
    for (auto& B: factor_starts[(A.second + 1)%n]) {
      int B_len = B.second - B.first + 1 + n * (B.second < B.first);
      // One empty 120-drome.
      if (A_len + B_len == n) {
        return true;
      }
      // No empty 120-drome factors.
      for (auto& C: factor_starts[(B.second+1)%n]) {
        int C_len = C.second - C.first + 1 + n * (C.second < C.first);
        if (A_len + B_len + C_len == n) {
          return true;
        }
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_case_8a_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<2*MAX_BND>& sixty_factors, const FactorArray<2*MAX_BND>& onetwenty_factors) {
  // At_60(A) Bt_120(B) C where C is palindrome.
  int n = P.size();
  vector<vector<Factor>> sixty_factor_starts(n);
  for (auto& f: sixty_factors) {
    sixty_factor_starts[f.first].push_back(f);
  }
  vector<vector<Factor>> onetwenty_factor_starts(n);
  for (auto& f: onetwenty_factors) {
    onetwenty_factor_starts[f.first].push_back(f);
  }

  // Factorizations with non-empty palindrome factor.
  for (auto& C: palin_factors) {
    int C_len = C.second - C.first + 1 + n * (C.second < C.first);
    // Empty 60-drome factor.
    for (auto& B: onetwenty_factor_starts[(C.second + 1)%n]) {
      int B_len = B.second - B.first + 1 + n * (B.second < B.first);
      if (B_len + C_len == n) {
        return true;
      }
    }
    // No empty 60-drome factor
    for (auto& A: sixty_factor_starts[(C.second + 1)%n]) {
      int A_len = A.second - A.first + 1 + n * (A.second < A.first);
      // Empty 120-drome factor.
      if (A_len + C_len == n) {
        return true;
      }
      // No empty 120-drome factor.
      for (auto& B: onetwenty_factor_starts[(A.second + 1)%n]) {
        int B_len = B.second - B.first + 1 + n * (B.second < B.first);
        if (A_len + B_len + C_len == n) {
          return true;
        }
      }
    }
  }
  // Factorizations with empty palindrome factor and non-empty 60-drome factor.
  for (auto& A: sixty_factors) {
    int A_len = A.second - A.first + 1 + n * (A.second < A.first);
    for (auto& B: sixty_factor_starts[(A.second + 1)%n]) {
      int B_len = B.second - B.first + 1 + n * (B.second < B.first);
      if (A_len + B_len == n) {
        return true;
      }
    }
  }
  return false;
}

bool IsohedralChecker::has_case_8b_tiling(const boundaryword& P, const FactorArray<2*MAX_BND>& palin_factors, const FactorArray<2*MAX_BND>& sixty_factors, const FactorArray<2*MAX_BND>& onetwenty_factors) {
  // At_60(A) B Ct_120(C)
  int n = P.size();
  vector<vector<Factor>> sixty_factor_starts(n);
  for (auto& f: sixty_factors) {
    sixty_factor_starts[f.first].push_back(f);
  }
  vector<vector<Factor>> onetwenty_factor_starts(n);
  for (auto& f: onetwenty_factors) {
    onetwenty_factor_starts[f.first].push_back(f);
  }

  // Factorizations with non-empty palindrome factor.
  for (auto& B: palin_factors) {
    int B_len = B.second - B.first + 1 + n * (B.second < B.first);
    // Empty 120-drome factor.
    for (auto& A: sixty_factor_starts[(B.second + 1)%n]) {
      int A_len = A.second - A.first + 1 + n * (A.second < A.first);
      if (A_len + B_len == n) {
        return true;
      }
    }
    // No empty 120-drome factor
    for (auto& C: onetwenty_factor_starts[(B.second + 1)%n]) {
      int C_len = C.second - C.first + 1 + n * (C.second < C.first);
      // Empty 60-drome factor.
      if (B_len + C_len == n) {
        return true;
      }
      // No empty 60-drome factor.
      for (auto& A: sixty_factor_starts[(C.second + 1)%n]) {
        int A_len = A.second - A.first + 1 + n * (A.second < A.first);
        if (A_len + B_len + C_len == n) {
          return true;
        }
      }
    }
  }
  // Factorizations with empty palindrome factor and no empty 120-drome factor.
  for (auto& C: onetwenty_factors) {
    int C_len = C.second - C.first + 1 + n * (C.second < C.first);
    for (auto& A: sixty_factor_starts[(C.second + 1)%n]) {
      int A_len = A.second - A.first + 1 + n * (A.second < A.first);
      if (A_len + C_len == n) {
        return true;
      }
    }
  }
  return false;

}


bool IsohedralChecker::has_isohedral_tiling(const boundaryword &P) {
  FactorPairArray<MAX_BND*MAX_BND> mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P);
  FactorArray<2*MAX_BND> palin_factors = admissible_rotadrome_factors(P, 180);
  if (has_half_turn_tiling(P, mirror_factor_pairs, palin_factors)) return true;

  FactorArray<2*MAX_BND> mirror_factors = admissible_mirror_factors(P);
  if (has_translation_tiling(P, mirror_factors)) return true;
  
  FactorArray<2*MAX_BND> ninety_factors = admissible_rotadrome_factors(P, 90);
  if (has_quarter_turn_tiling(P, ninety_factors, palin_factors)) {
    return true;
  }

  FactorArray<MAX_BND*MAX_BND> reflect_square_factors = admissible_reflect_square_factors(P);
  if (has_type_1_reflection_tiling(P, reflect_square_factors, mirror_factor_pairs)) {
    return true;
  }
  if (has_type_2_reflection_tiling(P, mirror_factors)) {
    return true;
  }

  if (has_type_1_half_turn_reflection_tiling(P, mirror_factor_pairs, palin_factors, reflect_square_factors)) {
    return true;
  }
  if (has_type_2_half_turn_reflection_tiling(P, palin_factors)) {
    return true;
  }

  FactorArray<2*MAX_BND> onetwenty_factors = admissible_rotadrome_factors(P, 120);
  if (has_case_7_tiling(P, onetwenty_factors)) {
    return true;
  }

  FactorArray<2*MAX_BND> sixty_factors = admissible_rotadrome_factors(P, 60);
  if (has_case_8a_tiling(P, palin_factors, sixty_factors, onetwenty_factors)) {
    return true;
  }
  if (has_case_8b_tiling(P, palin_factors, sixty_factors, onetwenty_factors)) {
    return true;
  }
  return false;
}

