#include "boundary.h"

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
  edges = getUniqueTileEdges(shape);
  cout << "New edges\n";
  for (auto& e: edges) {
    cout << e.first << "->" << e.second << "\n";
  }
}
