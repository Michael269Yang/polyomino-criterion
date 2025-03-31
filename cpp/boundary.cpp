#include "boundary.h"

boundaryword operator+(const boundaryword& lhs,
                       const boundaryword& rhs) {
  boundaryword result = lhs;
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}


boundaryword slice(const boundaryword& boundary, int start, int end) {
  if (end == -1) {
    end = boundary.size();
  }

  auto first = boundary.begin() + start;
  auto last = boundary.begin() + end;
  return boundaryword(first, last);
}
