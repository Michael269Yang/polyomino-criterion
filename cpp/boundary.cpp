#include <bitset>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <utility>

//#include <boost/functional/hash.hpp>

using std::cout;

enum GridType {
  NOGRID = -1,
  OMINO = 0,
};

template<typename coord>
class point
{
public:
	point() 
		: x_ {}, y_ {}
	{}
	template<typename ocoord>
	point( const point<ocoord>& other ) 
		: x_ { static_cast<coord>( other.x_ ) }
		, y_ { static_cast<coord>( other.y_ ) }
	{}
	point( coord x, coord y ) 
		: x_ { x }, y_ { y }
	{}

	coord getX() const { return x_; }
	coord getY() const { return y_; }

	template<typename ocoord>
	bool operator ==( const point<ocoord>& other ) const
	{
		return (x_ == other.x_) && (y_ == other.y_);
	}

	template<typename ocoord>
	bool operator !=( const point<ocoord>& other ) const
	{
		return (x_ != other.x_) || (y_ != other.y_);
	}

	template<typename ocoord>
	bool operator <( const point<ocoord>& other ) const
	{
		return (y_ < other.y_) || ((y_ == other.y_) && (x_ < other.x_));
	}

	template<typename ocoord>
	bool operator <=( const point<ocoord>& other ) const
	{
		return (y_ < other.y_) || ((y_ == other.y_) && (x_ <= other.x_));
	}

	template<typename ocoord>
	point<coord> operator +( const point<ocoord>& other ) const
	{
		return { coord( x_ + other.x_ ), coord( y_ + other.y_ ) };
	}

	template<typename ocoord>
	point<coord> operator -( const point<ocoord>& other ) const
	{
		return { coord( x_ - other.x_ ), coord( y_ - other.y_ ) };
	}

	point<coord> operator -() const
	{
		return { coord( -x_ ), coord( -y_ ) };
	}

	template<typename ocoord>
	point<coord>& operator =( const point<ocoord>& other )
	{
		x_ = other.x_;
		y_ = other.y_;
		return *this;
	}

	template<typename ocoord>
	point<coord>& operator +=( const point<ocoord>& other )
	{
		x_ += other.x_;
		y_ += other.y_;
		return *this;
	}

	template<typename ocoord>
	point<coord>& operator -=( const point<ocoord>& other ) 
	{
		x_ -= other.x_;
		y_ -= other.y_;
		return *this;
	}

	coord x_;
	coord y_;
};

template<typename coord>
inline size_t hash_value( const point<coord>& p )
{
	return p.hash();
}

template<typename coord>
inline std::ostream& operator <<( std::ostream& os, const point<coord>& p )
{
	return os << '<' << int(p.x_) << ',' << int(p.y_) << '>';
}

template<>
inline std::ostream& operator <<( std::ostream& os, const point<double>& p )
{
	return os << '<' << p.x_ << ',' << p.y_ << '>';
}

template<typename grid>
class Shape
{
public:
	using coord_t = typename grid::coord_t;
	using point_t = typename grid::point_t;
	using xform_t = typename grid::xform_t;

	Shape()
		: pts_ {}
	{}
	Shape( const Shape<grid>& other )
		: pts_ { other.pts_ }
	{}
	Shape( const Shape<grid>& other, const point_t& v )
	{
		for( auto p : other ) {
			pts_.push_back( p + v );
		}
	}

	size_t size() const
	{
		return pts_.size();
	}

	void add( const coord_t& x, const coord_t& y )
	{
		pts_.emplace_back( x, y );
	}
	void add( const point_t& p )
	{
		pts_.push_back( p );
	}
	void add( const Shape<grid>& other )
	{
		for( auto p : other ) {
			pts_.push_back( p );
		}
	}

	void complete()
	{
		pts_.sort();
	}

	Shape& operator =( const Shape& other ) 
	{
		if( pts_.size() == other.pts_.size() ) {
			std::copy( other.pts_.begin(), other.pts_.end(), pts_.begin() );
		} else {
			pts_.clear();
			std::copy( other.pts_.begin(), other.pts_.end(), 
				std::back_inserter( pts_ ) );
		}
		return *this;
	}

	void reset()
	{
		pts_.clear();
	}
	void reset( const Shape& other, const xform_t& T )
	{
		if( pts_.size() == other.pts_.size() ) {
			// If same size, overwrite without reallocating.
			auto i = pts_.begin();
			auto j = other.pts_.begin();
			while( i != pts_.end() ) {
				*i = T * (*j);
				++i;
				++j;
			}
		} else {
			pts_.clear();
			for( auto p : other ) {
				pts_.push_back( T * p );
			}
		}
		complete();
	}
	void translate( const point_t& dp )
	{
		// Translating should never affect the order.
		for( auto& p : pts_ ) {
			p += dp;
		}
	}
	// Reset with translation
	void reset( const Shape& other, const point_t& dp )
	{
		pts_.clear();
		for( auto p : other ) {
			pts_.push_back( p + dp );
		}
	}

	bool intersects( const Shape<grid>& other ) const;
	bool operator ==( const Shape<grid>& other ) const
	{
		return std::equal(
			pts_.begin(), pts_.end(), other.pts_.begin(), other.pts_.end() );
	}
	bool operator !=( const Shape<grid>& other ) const
	{
		return (*this) != other;
	}

	void getSymmetries( std::vector<xform_t>& syms ) const
	{
		Shape<grid> other;
		syms.clear();
		for( size_t idx = 1; idx < grid::num_orientations; ++idx ) {
			other.reset( *this, grid::orientations[idx] );
			if( equivalent( other ) ) {
				syms.push_back( grid::orientations[idx] );
			}
		}
	}

	// Are these shapes equivalent under translation?
	bool equivalent( const Shape<grid>& other ) const
	{
		if( size() != other.size() ) {
			return false;
		}

		if( !grid::translatable( pts_.front(), other.pts_.front() ) ) {
			return false;
		}
		
		point_t d = pts_.front() - other.pts_.front();

		auto i = pts_.begin();
		auto j = other.pts_.begin();

		while( i != pts_.end() ) {
			if( (*i) != ((*j)+d) ) {
				return false;
			}
			++i;
			++j;
		}

		return true;
	}

	// Move this shape so its minimum point lies at an origin of the grid.
	void untranslate()
	{
		point_t p = pts_.front();
		point_t v = grid::getOrigin( p ) - p;

		for( auto& sp : pts_ ) {
			sp += v;
		}
	}

	int compare( const Shape<grid>& other ) const
	{
		auto i = pts_.begin();
		auto j = other.pts_.begin();

		while( true ) {
			bool ei = (i == pts_.end());
			bool ej = (j == other.pts_.end());

			if( ei && ej ) {
				return 0;
			} else if( ei ) {
				return -1;
			} else if( ej ) {
				return 1;
			} else if( *i < *j ) {
				return -1;
			} else if( *j < *i ) {
				return 1;
			}

			++i;
			++j;
		}
	}
	bool operator <( const Shape<grid>& other ) const
	{
		return compare( other ) < 0;
	}

	void getHaloAndBorder( Shape<grid>& halo, Shape<grid>& border ) const;
	void getEdgeHalo( Shape<grid>& halo ) const;
	bool simplyConnected() const;

	auto begin()
	{ return pts_.begin(); }
	auto end()
	{ return pts_.end(); }

	auto begin() const
	{ return pts_.begin(); }
	auto end() const
	{ return pts_.end(); }

	void debug() const;

private:
	std::list<point_t> pts_;
};

template<typename coord>
class xform
{
public:
	xform()
		: a_ { 1 } , b_ { 0 } , c_ { 0 }
		, d_ { 0 } , e_ { 1 } , f_ { 0 }
	{}
	xform( coord a, coord b, coord c, coord d, coord e, coord f )
		: a_ { a } , b_ { b } , c_ { c }
		, d_ { d } , e_ { e } , f_ { f }
	{}
	template<typename ocoord>
	xform( const xform<ocoord>& other )
		: a_ { static_cast<coord>( other.a_ ) } 
		, b_ { static_cast<coord>( other.b_ ) } 
		, c_ { static_cast<coord>( other.c_ ) }
		, d_ { static_cast<coord>( other.d_ ) } 
		, e_ { static_cast<coord>( other.e_ ) } 
		, f_ { static_cast<coord>( other.f_ ) }
	{}

	template<typename ocoord>
	xform<coord>& operator =( const xform<ocoord>& other )
	{
		a_ = other.a_;
		b_ = other.b_;
		c_ = other.c_;
		d_ = other.d_;
		e_ = other.e_;
		f_ = other.f_;
		return *this;
	}

	template<typename ocoord>
	point<coord> operator *( const point<ocoord>& p ) const
	{
		return { coord( a_ * p.x_ + b_ * p.y_ + c_ ),
				 coord( d_ * p.x_ + e_ * p.y_ + f_ ) };
	}
	template<typename ocoord>
	xform<coord> operator *( const xform<ocoord>& other ) const
	{
		return {
			coord( a_*other.a_ + b_*other.d_ ),
			coord( a_*other.b_ + b_*other.e_ ),
			coord( a_*other.c_ + b_*other.f_ + c_ ),

			coord( d_*other.a_ + e_*other.d_ ),
			coord( d_*other.b_ + e_*other.e_ ),
			coord( d_*other.c_ + e_*other.f_ + f_ ) };
	}
	template<typename ocoord>
	xform<coord> translate( const point<ocoord>& pt ) const
	{
		return { a_, b_, coord( c_ + pt.x_ ), d_, e_, coord( f_ + pt.y_ ) };
	}

	xform<coord> invert() const
	{
		// Must be +-1, no need to take reciprocal
		coord det( a_ * e_ - b_ * d_ );

		return { 
			coord( e_ * det ), 
			coord( -b_ * det ),
			coord( (b_ * f_ - c_ * e_) * det ),

			coord( -d_ * det ), 
			coord( a_ * det ),
			coord( (c_ * d_ - a_ * f_) * det ) };
	}

	template<typename ocoord>
	xform<coord>& operator +=( const point<ocoord>& pt ) 
	{
		c_ += pt.x_;
		f_ += pt.y_;
		return *this;
	}

	template<typename ocoord>
	bool operator ==( const xform<ocoord>& other ) const
	{
		return (a_ == other.a_) 
			&& (b_ == other.b_) 
			&& (c_ == other.c_) 
			&& (d_ == other.d_) 
			&& (e_ == other.e_) 
			&& (f_ == other.f_);
	}
	template<typename ocoord>
	bool operator !=( const xform<ocoord>& other ) const
	{
		return !(*this == other);
	}

	bool isIdentity() const
	{
		return (a_==1) && (b_==0) && (c_==0) && (d_==0) && (e_==1) && (f_==0);
	}
	bool isTranslation() const
	{
		return (a_==1) && (b_==0) && (d_==0) && (e_==1);
	}
	bool isHalfturn() const
	{
		return (a_==-1) && (b_==0) && (d_==0) && (e_==-1);
	}
	coord det() const
	{
		return a_*e_ - b_*d_;
	}

	/*size_t hash() const
	{
		size_t res = 0;
		boost::hash_combine( res, a_ );
		boost::hash_combine( res, b_ );
		boost::hash_combine( res, c_ );
		boost::hash_combine( res, d_ );
		boost::hash_combine( res, e_ );
		boost::hash_combine( res, f_ );
		return res;
	}*/

	coord a_;
	coord b_;
	coord c_;
	coord d_;
	coord e_;
	coord f_;
};

template<typename coord>
class OminoGrid
{
public:
	using coord_t = coord;
	using point_t = point<coord>;
	using xform_t = xform<coord>;
	using edge_t = std::pair<point_t, point_t>;

    enum TileType {
		INVALID = -1,
		SQUARE = 0
    };

    enum TileShape {
		SQUARE_SHAPE = 0
	};

public:
	inline static GridType grid_type = OMINO;

	// Number of transivity classes of tiles under translation
    inline static size_t num_tile_types = 1; 
	// Number of distinct shapes 
    inline static size_t num_tile_shapes = 1;
	// What tile type is the tile indexed by p?
	inline static TileType getTileType( const point_t& p )
	{
		return SQUARE;
	}
	// What shape tile is at this position?
	inline static TileShape getTileShape( const point_t& p )
	{
		return SQUARE_SHAPE;
	}
	// Get the origin point 
	inline static point_t getOrigin( const point_t& p )
	{
		return { 0, 0 };
	}

	inline static size_t numNeighbours( const point_t& p )
	{
		return 8;
	}

	static const point<int8_t> *getNeighbourVectors( const point_t& p )
	{
		return all_neighbours;
	}

	static size_t numEdgeNeighbours( const point_t& p )
	{
		return 4;
	}

	static const point<int8_t> *getEdgeNeighbourVectors( const point_t& p )
	{
		return edge_neighbours;
	}

	static bool translatable( const point_t& p, const point_t& q )
	{
		return true;
	}

	// Functions to assist with rendering

	// Get points with integer coordinates corresponding to the *vertices*
	// of the cell indexed by p.  These can really be in any coordinate
	// system whatsoever -- they just need to be in one-to-one correspondence
	// with the actual vertices, so that they can be compared exactly.
	static std::vector<point_t> getCellVertices( const point_t& p )
	{
		return {
			p, 
			p + point_t { 1, 0 },
			p + point_t { 1, 1 },
			p + point_t { 0, 1 } };
	}

	// Convert a vertex as given by the previous function into a 2D
	// point that's compatible with the matrices giving the symmetries
	// of the grid.
	static point<double> vertexToGrid( const point_t& pt ) 
	{
		return point<double>( 
			static_cast<double>(pt.getX()) - 0.5,
			static_cast<double>(pt.getY()) - 0.5 );
	}

	// Any final cleanup of a point given by the previous function.
	static point<double> gridToPage( const point<double>& pt )
	{
		return pt;
	}

	static const point_t origins[1];

	static const size_t num_orientations;
	static const xform<int8_t> orientations[8];
	
	static const point<int8_t> all_neighbours[8];
	static const point<int8_t> edge_neighbours[4];
};

template<typename coord>
using edge = std::pair<point<coord>, point<coord>>;
template<typename coord>
using edgeset = std::set<edge<coord>>;

template<typename coord>
void pe(const edge<coord>& e) {
  cout << e.first << "->" << e.second << "\n";
}

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

    if (pt.getX() < bottomLeft.getX() || 
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

int main() {
  cout << "Hello World\n";

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
