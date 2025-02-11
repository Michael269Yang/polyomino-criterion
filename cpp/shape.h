#ifndef SHAPE_H
#define SHAPE_H

#include <list>

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

#endif // SHAPE_H

