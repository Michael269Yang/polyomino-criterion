#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "geom.h"
#include "hexgrid.h"
#include "iamondgrid.h"
#include "ominogrid.h"
#include "shape.h"

#include <map>
#include <set>
#include <utility>

using boundaryword = std::vector<std::pair<int, int>>;

boundaryword operator+(const boundaryword& lhs,
                       const boundaryword& rhs);
boundaryword slice(const boundaryword& boundary, int start, int end=-1);

template<typename coord>
using edge = std::pair<point<coord>, point<coord>>;
template<typename coord>
using edgeset = std::set<edge<coord>>;
template<typename coord>
using edgemap = std::map<point<coord>, std::set<point<coord>>>;

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

// This code currently only works for the hex grid.
template <typename coord>
boundaryword getBoundaryWord(const Shape<HexGrid<coord>>& shape) {
  using edge_t = edge<coord>;
  using edgeset_t = edgeset<coord>;
  using edgemap_t = edgemap<coord>;
  using point_t = typename HexGrid<coord>::point_t;

  edgemap_t neighbours;
  point_t bottomLeft(300, 300);
  for (auto& edge: getUniqueTileEdges(shape)) {
    neighbours[edge.first].insert(edge.second);
    neighbours[edge.second].insert(edge.first);

    bottomLeft = std::min(bottomLeft, edge.first);
    bottomLeft = std::min(bottomLeft, edge.second);
  }

  point_t start = bottomLeft;
  point_t cur = start;

  std::pair<int, int> NW = {-2, 1};
  std::pair<int, int> NE = {1, 1};
  std::pair<int, int> SW = {-1, -1};
  std::pair<int, int> SE = {2, -1};
  std::pair<int, int> U = {-1, 2};

  boundaryword boundary;

  // From the bottom left we can only go NW or U. Prefer NW.
  point_t nw = start + point_t(-2, 1);
  if (neighbours[start].find(nw) != neighbours[start].end()) {
    cur = nw;
    neighbours[nw].erase(start);
    boundary.push_back(NW);
  } else {
    cur = start + point_t(-1, 2);
    neighbours[cur].erase(start);
    boundary.push_back({-1, 2});
  }

  /*std::map<point_t, char> edgeToLetter = {
    {point_t(-1, 2), 'U'}, {point_t(1, -2), 'D'}, {point_t(1, 1), 'R'},
    {point_t(2, -1), 'r'}, {point_t(-2, 1), 'L'}, {point_t(-1, -1), 'l'}
  };*/
  while (cur != start) {
    point_t next = *neighbours[cur].begin();
    point_t edgeDir = next - cur;
    neighbours[next].erase(cur);
    cur = next;
    boundary.push_back({edgeDir.getX(), edgeDir.getY()});
  }
  return boundary;
}

// This code currently only works for the iamond grid.
template <typename coord>
boundaryword getBoundaryWord(const Shape<IamondGrid<coord>>& shape) {
  using edge_t = edge<coord>;
  using edgeset_t = edgeset<coord>;
  using edgemap_t = edgemap<coord>;
  using point_t = typename HexGrid<coord>::point_t;

  edgemap_t neighbours;
  point_t bottomLeft(300, 300);
  for (auto& edge: getUniqueTileEdges(shape)) {
    neighbours[edge.first].insert(edge.second);
    neighbours[edge.second].insert(edge.first);

    bottomLeft = std::min(bottomLeft, edge.first);
    bottomLeft = std::min(bottomLeft, edge.second);
  }

  point_t start = bottomLeft;
  point_t cur = start;

  std::pair<int, int> E = {3, 0};
  std::pair<int, int> NE = {0, 3};
  std::pair<int, int> NW = {-3, 3};
  std::pair<int, int> W = {-3, 0};
  std::pair<int, int> SW = {0, -3};
  std::pair<int, int> SE = {3, -3};

  boundaryword boundary;

  // From the bottom left we can only go NW or NE. Prefer NW.
  point_t nw = start + point_t(NW.first, NW.second);
  if (neighbours[start].find(nw) != neighbours[start].end()) {
    cur = nw;
    neighbours[nw].erase(start);
    boundary.push_back(NW);
  } else {
    cur = start + point_t(NE.first, NE.second);
    neighbours[cur].erase(start);
    boundary.push_back(NE);
  }

  while (cur != start) {
    point_t next = *neighbours[cur].begin();
    point_t edgeDir = next - cur;
    neighbours[next].erase(cur);
    cur = next;
    boundary.push_back({edgeDir.getX(), edgeDir.getY()});
  }
  return boundary;
}


template <typename coord>
boundaryword getBoundaryWord(const Shape<OminoGrid<coord>>& shape) {
  using point_t = typename OminoGrid<coord>::point_t;

  const int E = 0;
  const int N = 1;
  const int W = 2;
  const int S = 3;

  const std::vector<std::pair<int, int>> dirs = 
    {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

  boundaryword ret;

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

  // Traverse in clockwise direction. 
  // Since we're at bottom left point, next point should be up. 
  // Start at next point and keep traversing boundary until reach start aka bottomLeft.
  point_t cur(bottomLeft.getX(), bottomLeft.getY() + 1);
  int8_t prevDir = N;
  ret.push_back({0, 1});
  while (cur != bottomLeft) {
    // If we're currently going straight, always prioritize going left, then 
    // continuing straight, then going right.
    int8_t rightDir = (prevDir + 3) % 4;
    int8_t leftDir = (prevDir + 1) % 4;

    int8_t prevDirX = (prevDir == E || prevDir == W) ? (prevDir == E ? 1 : -1) : 0;
    int8_t prevDirY = (prevDir == N || prevDir == S) ? (prevDir == N ? 1 : -1) : 0;
    point_t next;
    if (vertexToEdges[cur].test(leftDir)) {
      // Turn left.
      // 90 degree rotation of direction (x, y) is (-y, x).
      prevDir = leftDir;
      next = cur + point_t(-prevDirY, prevDirX);
      ret.push_back(dirs[leftDir]);
    } else if (vertexToEdges[cur].test(prevDir)) {
      // Go straight. 
      // 0 degree rotation of direction (x, y) is (x, y).
      next = cur + point_t(prevDirX, prevDirY);
      ret.push_back(dirs[prevDir]);
    } else {
      // Turn right. 
      // 270 degree rotation of direction (x, y) is (y, -x).
      prevDir = rightDir;
      next = cur + point_t(prevDirY, -prevDirX);
      ret.push_back(dirs[rightDir]);
    }

    cur = next;
  }
  return ret;
}

#endif // BOUNDARY_H

