#include "boundary.h"
#include "ominogrid.h"

#include <iostream>

using namespace std;

int main() {
  Shape<OminoGrid<int>> shape;
  shape.add(0, 0);
  auto edges = getUniqueTileEdges(shape);
  if (edges.empty()) {
    cout << "Edges is empty\n";
  }
  else for (auto& e: edges ) {
    cout << "(" << e.first.getX() << ", " << e.first.getY() << ") to (" << e.second.getX() << ", " << e.second.getY() << ")\n";
  }

  shape = Shape<OminoGrid<int>>();
  shape.add(0, 0);
  shape.add(1, 0);
  shape.add(2, 0);
  shape.add(2, -1);
  shape.add(3, 0);
  shape.add(2, 1);
  cout << getBoundaryWord(shape) << "\n";

  cout << "Testing for square\n";
  shape = Shape<OminoGrid<int>>();
  shape.add(0, 0);
  cout << getBoundaryWord(shape) << "\n";

  cout << "Testing 5-omino\n";
  shape = Shape<OminoGrid<int>>();
  shape.add(0, 0);
  shape.add(0, 1);
  shape.add(1, 1);
  shape.add(2, 1);
  shape.add(2, 0);
  cout << "Boundary: " << getBoundaryWord(shape) << "\n";
}
