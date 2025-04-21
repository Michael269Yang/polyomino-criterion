#include "boundary.h"
#include "hexgrid.h"
#include "iamondgrid.h"
#include "ominogrid.h"

#include <cassert>
#include <iostream>
#include <map>
#include <string>

using namespace std;

int main() {
  map<pair<int, int>, char> ominoMap{
    {{1, 0}, 'E'},
    {{-1, 0}, 'W'},
    {{0, 1}, 'N'},
    {{0, -1}, 'S'}
  };
  Shape<OminoGrid<int>> ominoShape;
  ominoShape.add(0, 0);
  ominoShape.add(0, 1);
  ominoShape.add(1, 1);
  boundaryword ominoBoundary = getBoundaryWord(ominoShape);
  string got;
  for (auto& p: ominoBoundary) {
    got += ominoMap[p];
  }
  string want = "NNEESWSW";
  if (got == want) {
    cout << "Omino boundary test case passed.\n";
  } else {
    cout << "Omino boundary test failed.\n";
    cout << "Got: " << got << " Want: " << want << "\n";
  }

  map<pair<int, int>, char> hexMap{
    {{-1, 2}, 'U'},
    {{1, -2}, 'D'},
    {{1, 1}, 'R'},
    {{2, -1}, 'r'},
    {{-2, 1}, 'L'},
    {{-1, -1}, 'l'}
  };
  Shape<HexGrid<int>> hexShape;
  hexShape.add(0, 0);
  hexShape.add(0, 1);
  hexShape.add(1, 0);
  boundaryword hexBoundary = getBoundaryWord(hexShape);
  got = "";
  for (auto& p: hexBoundary) {
    got += hexMap[p];
  }
  want = "LURURrDrDlLl";
  if (got == want) {
    cout << "Hex boundary test case passed.\n";
  } else {
    cout << "Hex boundary test failed.\n";
    cout << "Got: " << got << " Want: " << want << "\n";
  }

  map<pair<int, int>, string> iamondMap{
    {{3, 0}, "E"},
    {{0, 3}, "nE"},
    {{-3, 3}, "nW"},
    {{-3, 0}, "W"},
    {{0, -3}, "sW"},
    {{3, -3}, "sE"}
  };
  Shape<IamondGrid<int>> iamondShape;
  iamondShape.add(0, 0);
  iamondShape.add(1, 1);
  boundaryword iamondBoundary = getBoundaryWord(iamondShape);
  got = "";
  for (auto& p: iamondBoundary) {
    got += iamondMap[p];
  }
  want = "nEEsWW";
  if (got == want) {
    cout << "Iamond boundary test case passed.\n";
  } else {
    cout << "Iamond boundary test failed.\n";
    cout << "Got: " << got << " Want: " << want << "\n";
  }
}
