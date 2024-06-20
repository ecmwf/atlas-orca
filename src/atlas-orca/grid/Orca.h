/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"
#include "atlas-orca/util/PointIJ.h"

namespace atlas {
class Mesh;
namespace orca {
class OrcaPeriodicity;
}
}  // namespace atlas
namespace eckit {
class PathName;
}


namespace atlas::grid::detail::grid {

class Orca final : public Grid {
private:
    struct ComputePointXY {
        explicit ComputePointXY( const Orca& grid ) : grid_( grid ) {}

        void operator()( idx_t i, idx_t j, PointXY& point ) { grid_.xy( i, j, point.data() ); }

        const Orca& grid_;
    };

    struct ComputePointLonLat {
        explicit ComputePointLonLat( const Orca& grid ) : grid_( grid ) {}

        void operator()( idx_t i, idx_t j, PointLonLat& point ) { grid_.lonlat( i, j, point.data() ); }

        const Orca& grid_;
    };

    template <typename Base, typename ComputePoint>
    class OrcaIterator : public Base {
        const Orca& grid_;
        idx_t ibegin_;
        idx_t iend_;
        idx_t jbegin_;
        idx_t jend_;
        idx_t i_;
        idx_t j_;
        typename Base::value_type point_;
        ComputePoint compute_point_;

    public:
        explicit OrcaIterator( const Orca& grid, bool begin = true ) :
            grid_( grid ),
            ibegin_( -grid.haloWest() ),
            iend_( grid.nx() + grid.haloEast() ),
            jbegin_( -grid.haloSouth() ),
            jend_( grid.ny() + grid.haloNorth() ),
            i_( -grid.haloWest() ),
            j_( -grid.haloSouth() ),
            compute_point_{ grid_ } {
            if ( not begin ) {
                i_ = ibegin_;
                j_ = jend_;
            }
            if ( j_ != jend_ && i_ < iend_ ) {
                compute_point_( i_, j_, point_ );
            }
        }

        bool next( typename Base::value_type& point ) override {
            if ( j_ < jend_ && i_ < iend_ ) {
                compute_point_( i_++, j_, point );
                if ( i_ == iend_ ) {
                    ++j_;
                    i_ = ibegin_;
                }
                return true;
            }
            return false;
        }

        const typename Base::reference operator*() const override { return point_; }

        const Base& operator++() override {
            ++i_;
            if ( i_ == iend_ ) {
                ++j_;
                i_ = ibegin_;
            }
            if ( j_ != jend_ && i_ < iend_ ) {
                compute_point_( i_, j_, point_ );
            }
            return *this;
        }


        const Base& operator+=( typename Base::difference_type distance ) override {
            idx_t d = distance;
            while ( j_ != jend_ && d >= ( iend_ - i_ ) ) {
                d -= ( iend_ - i_ );
                ++j_;
                i_ = ibegin_;
            }
            i_ += d;
            if ( j_ != jend_ && i_ < iend_ ) {
                compute_point_( i_, j_, point_ );
            }
            return *this;
        }

        typename Base::difference_type distance( const Base& other ) const override {
            const auto& _other               = static_cast<const OrcaIterator&>( other );
            typename Base::difference_type d = 0;
            idx_t j                          = j_;
            idx_t i                          = i_;
            while ( j < _other.j_ ) {
                d += iend_ - i;
                ++j;
                i = 0;
            }
            d += _other.i_;
            return d;
        }

        bool operator==( const Base& other ) const override {
            return j_ == static_cast<const OrcaIterator&>( other ).j_ &&
                   i_ == static_cast<const OrcaIterator&>( other ).i_;
        }

        bool operator!=( const Base& other ) const override { return !( *this == other ); }

        std::unique_ptr<Base> clone() const override {
            auto result    = new OrcaIterator( grid_, false );
            result->i_     = i_;
            result->j_     = j_;
            result->point_ = point_;
            return std::unique_ptr<Base>( result );
        }
    };

public:
    using IteratorXY     = OrcaIterator<Grid::IteratorXY, ComputePointXY>;
    using IteratorLonLat = OrcaIterator<Grid::IteratorLonLat, ComputePointLonLat>;

public:  // methods
    static std::string static_type();

    /// Constructor taking a configuration (spec)
    explicit Orca( const Config& );

    /// Constructor taking a name/uid and a configuration (spec)
    Orca( const std::string& name_or_uid, const Config& );

    idx_t size() const override;

    Spec spec() const override;

    idx_t nx() const { return nx_; }
    idx_t ny() const { return ny_; }

    std::string name() const override;
    std::string type() const override;

    const PointXY& xy( idx_t i, idx_t j ) const { return points_[index( i, j )]; }

    void xy( idx_t i, idx_t j, double crd[] ) const {
        const PointXY& p = xy( i, j );
        crd[0]           = p[0];
        crd[1]           = p[1];
    }

    PointLonLat lonlat( idx_t i, idx_t j ) const { return xy( i, j ); }

    ATLAS_ALWAYS_INLINE gidx_t index( idx_t i, idx_t j ) const { return ( imin_ + i ) + ( jmin_ + j ) * jstride_; }

    void lonlat( idx_t i, idx_t j, double crd[] ) const { xy( i, j, crd ); }

    bool water( idx_t i, idx_t j ) const { return water_[( imin_ + i ) + ( jmin_ + j ) * jstride_]; }
    bool land( idx_t i, idx_t j ) const { return not water_[( imin_ + i ) + ( jmin_ + j ) * jstride_]; }
    bool ghost( idx_t i, idx_t j ) const { return ghost_[( imin_ + i ) + ( jmin_ + j ) * jstride_]; }
    bool invalidElement( idx_t i, idx_t j ) const { return invalid_element_[( imin_ + i ) + ( jmin_ + j ) * jstride_]; }

    gidx_t periodicIndex( idx_t i, idx_t j ) const;
    atlas::orca::PointIJ periodicIJ( idx_t i, idx_t j ) const;

    void index2ij( gidx_t gidx, idx_t& i, idx_t& j ) const {
        //gidx = jstride_ * (jmin_+j) + (imin_+i);
        j = static_cast<idx_t>( gidx / jstride_ - jmin_ );
        i = static_cast<idx_t>( gidx - ( jstride_ * ( j + jmin_ ) ) - imin_ );
    }


    std::unique_ptr<Grid::IteratorXY> xy_begin() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this ) );
    }
    std::unique_ptr<Grid::IteratorXY> xy_end() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this, false ) );
    }
    std::unique_ptr<Grid::IteratorLonLat> lonlat_begin() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this ) );
    }
    std::unique_ptr<Grid::IteratorLonLat> lonlat_end() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this, false ) );
    }

    int haloSouth() const { return halo_south_; }
    int haloNorth() const { return halo_north_; }
    int haloWest() const { return halo_west_; }
    int haloEast() const { return halo_east_; }

    size_t footprint() const override;

    Config meshgenerator() const override;

    Config partitioner() const override;

private:  // methods
    void print( std::ostream& ) const override;

    void hash( eckit::Hash& ) const override;

    std::string uid() const override;

    RectangularLonLatDomain lonlatBoundingBox() const override;

private:
    /// Grid name
    const std::string name_;

    /// Grid size
    idx_t nx_;
    idx_t ny_;
    idx_t halo_east_;
    idx_t halo_west_;
    idx_t halo_south_;
    idx_t halo_north_;
    idx_t nx_halo_;
    idx_t ny_halo_;

    /// Indices
    idx_t imin_;
    idx_t jmin_;
    idx_t istride_;
    idx_t jstride_;

    /// Storage of coordinate points
    std::vector<PointXY> points_;
    std::vector<bool> water_;
    std::vector<bool> ghost_;
    std::vector<bool> invalid_element_;

    /// Grid spec
    Spec spec_;

    std::unique_ptr<orca::OrcaPeriodicity> periodicity_;
};


}  // namespace atlas::grid::detail::grid
