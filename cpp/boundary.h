#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "geom.h"
#include "hexgrid.h"
#include "iamondgrid.h"
#include "kitegrid.h"
#include "ominogrid.h"
#include "shape.h"

#include <map>
#include <set>
#include <utility>

using boundaryword = std::vector<std::pair<int, int>>;

// Returns concatenation of 2 boundary words.
boundaryword operator+(const boundaryword& lhs,
                       const boundaryword& rhs);
// Returns slice of 0-indexed boundary [start, end). Returns start to end of the boundary word if end=-1.
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

// Given area based representation of a shape in the hex grid, returns
// a clockwise boundary word starting at the bottom left point.
template <typename coord>
boundaryword getBoundaryWord(const Shape<HexGrid<coord>>& shape) {
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
  std::pair<int, int> N = {-1, 2};

  boundaryword boundary;

  // From the bottom left we can only go NW or U. Prefer NW.
  point_t nw = start + point_t(NW.first, NW.second);
  if (neighbours[start].find(nw) != neighbours[start].end()) {
    cur = nw;
    neighbours[nw].erase(start);
    boundary.push_back(NW);
  } else {
    cur = start + point_t(N.first, N.second);
    neighbours[cur].erase(start);
    boundary.push_back(N);
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

// Given area based representation of a shape in the kite grid, returns
// a clockwise boundary word starting at the bottom left point.
template <typename coord>
boundaryword getBoundaryWord(const Shape<KiteGrid<coord>>& shape) {
  using edgemap_t = edgemap<coord>;
  using point_t = typename KiteGrid<coord>::point_t;

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
  std::pair<int, int> sNW = {-1, 1};
  std::pair<int, int> N = {-1, 2};
  std::pair<int, int> sNE = {0, 1};
  std::pair<int, int> NE = {1, 1};

  boundaryword boundary;

  // From the bottom left we pick greedily in the order NW, sNW, N, sNE, NE.
  point_t nw = start + point_t(NW.first, NW.second);
  point_t snw = start + point_t(sNW.first, sNW.second);
  point_t n = start + point_t(N.first, N.second);
  point_t sne = start + point_t(sNE.first, sNE.second);
  point_t ne = start + point_t(NE.first, NE.second);
  if (neighbours[start].find(nw) != neighbours[start].end()) {
    cur = nw;
    boundary.push_back(NW);
  } else if (neighbours[start].find(snw) != neighbours[start].end()) {
    cur = snw;
    boundary.push_back(sNW);
  } else if (neighbours[start].find(n) != neighbours[start].end()) {
    cur = n; 
    boundary.push_back(N);
  } else if (neighbours[start].find(sne) != neighbours[start].end()) {
    cur = sne;
    boundary.push_back(sNE);
  } else {
    cur = ne;
    boundary.push_back(NE);
  }
  neighbours[cur].erase(start);

  while (cur != start) {
    point_t next = *neighbours[cur].begin();
    point_t edgeDir = next - cur;
    neighbours[next].erase(cur);
    cur = next;
    boundary.push_back({edgeDir.getX(), edgeDir.getY()});
  }
  return boundary;
}



// Given area based representation of a shape in the iamond grid, returns
// a clockwise boundary word starting at the bottom left point.
template <typename coord>
boundaryword getBoundaryWord(const Shape<IamondGrid<coord>>& shape) {
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

  std::pair<int, int> NE = {0, 3};
  std::pair<int, int> NW = {-3, 3};

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


// Given area based representation of a shape in the omino grid, returns
// a clockwise boundary word starting at the bottom left point.
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

