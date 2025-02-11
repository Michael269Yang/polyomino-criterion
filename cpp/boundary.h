#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "geom.h"
#include "ominogrid.h"
#include "shape.h"

#include <map>
#include <set>
#include <utility>

template<typename coord>
using edge = std::pair<point<coord>, point<coord>>;
template<typename coord>
using edgeset = std::set<edge<coord>>;

template<typename grid>
edgeset<typename grid::coord_t> getUniqueTileEdges(const Shape<grid>& shape) {
  using coord_t = typename grid::coord_t;
  using edge_t = edge<coord_t>;
  using edgeset_t = edgeset<coord_t>;
  using point_t = typename grid::point_t;

  edgeset_t ret;
  for (const point_t& pt: shape) {
    std::vector<point_t> verts = grid::getCellVertices( pt );
    point_t last = verts.back();
    for (const point_t& cur: verts) {
      edge_t e { last, cur };
      edge_t opp { cur, last };
      if ( ret.find( opp ) == ret.end() ) {
        ret.insert( e );
      } else {
        ret.erase( opp );
      }
      last = cur;
    }
  }

  return ret;
}

template <typename coord>
edgeset<coord> getUniqueTileEdges(const Shape<OminoGrid<coord>>& shape) {
  using point_t = typename OminoGrid<coord>::point_t;

  const int E = 0;
  const int N = 1;
  const int W = 2;
  const int S = 3;

  edgeset<coord> ret;

  // bit 0 = N edge, bit 1 = E edge, bit 2 = S edge, bit 3 = W edge.
  std::map<point_t, std::bitset<4>> vertexToEdges;
  point_t bottomLeft(100, 100);

  for (const point_t& pt: shape) {
    point_t right = pt + point_t(1, 0);
    point_t up = pt + point_t(0, 1);
    point_t rightUp = pt + point_t(1, 1);

    // Has N and E edge.
    vertexToEdges[pt].set(N);
    vertexToEdges[pt].set(E);

    // Has N and W edge.
    vertexToEdges[right].set(N);
    vertexToEdges[right].set(W);

    // Has E and S edge.
    vertexToEdges[up].set(E);
    vertexToEdges[up].set(S);

    // Has S and W edge.
    vertexToEdges[rightUp].set(S);
    vertexToEdges[rightUp].set(W);

    if ((pt.getX() < bottomLeft.getX()) || 
      (pt.getX() == bottomLeft.getX() && pt.getY() < bottomLeft.getY())) {
      bottomLeft = pt;
    }
  }

  // Traverse in counter-clockwise direction. 
  // Since we're at bottom left point, next point should be to right. 
  // Start at next point and keep traversing boundary until reach start aka bottomLeft.
  point_t cur(bottomLeft.getX() + 1, bottomLeft.getY());
  int8_t prevDir = E;
  ret.insert({bottomLeft, cur});
  while (cur != bottomLeft) {
    // If we're currently going straight, always prioritize going right, then 
    // continuing straight, then going left.
    int8_t rightDir = (prevDir + 3) % 4;
    int8_t leftDir = (prevDir + 1) % 4;

    int8_t prevDirX = (prevDir == E || prevDir == W) ? (prevDir == E ? 1 : -1) : 0;
    int8_t prevDirY = (prevDir == N || prevDir == S) ? (prevDir == N ? 1 : -1) : 0;
    point_t next;
    if (vertexToEdges[cur].test(rightDir)) {
      // Turn right. 
      // 270 degree rotation of direction (x, y) is (y, -x).
      prevDir = rightDir;
      next = cur + point_t(prevDirY, -prevDirX);
    } else if (vertexToEdges[cur].test(prevDir)) {
      // Go straight. 
      // 0 degree rotation of direction (x, y) is (x, y).
      next = cur + point_t(prevDirX, prevDirY);
    } else {
      // Turn left.
      // 90 degree rotation of direction (x, y) is (-y, x).
      prevDir = leftDir;
      next = cur + point_t(-prevDirY, prevDirX);
    }

    ret.insert({cur, next});
    cur = next;
  }
  return ret;
}

#endif // BOUNDARY_H

