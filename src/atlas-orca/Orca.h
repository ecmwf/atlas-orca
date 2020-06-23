/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file Orca.h
/// @author Domokos Sarmany
/// @date January 2020

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"

namespace atlas {
class Mesh;
}

namespace eckit {
class PathName;
}

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class Orca final : public Grid {
private:
    struct ComputePointXY {
        ComputePointXY( const Orca& grid ) : grid_( grid ) {}

        void operator()( idx_t i, idx_t j, PointXY& point ) { grid_.xy( i, j, point.data() ); }

        const Orca& grid_;
    };

    struct ComputePointLonLat {
        ComputePointLonLat( const Orca& grid ) : grid_( grid ) {}

        void operator()( idx_t i, idx_t j, PointLonLat& point ) { grid_.lonlat( i, j, point.data() ); }

        const Orca& grid_;
    };

    template <typename Base, typename ComputePoint>
    class OrcaIterator : public Base {
    protected:
        const Orca& grid_;
        idx_t size_;
        idx_t i_;
        idx_t j_;
        typename Base::value_type point_;
        ComputePoint compute_point_;

    public:
        OrcaIterator( const Orca& grid, bool begin = true ) :
            grid_( grid ),
            size_( static_cast<idx_t>( grid_.points_.size() ) ),
            i_( 0 ),
            j_( 0 ),
            compute_point_{grid_} {
            compute_point_( i_, j_, point_ );
        }

        bool next( typename Base::value_type& point ) override {
            if ( j_ < grid_.ny() && i_ < grid_.nx() ) {
                compute_point_( i_++, j_, point );

                if ( i_ == grid_.nx() ) {
                    ++j_;
                    i_ = 0;
                }
                return true;
            }
            return false;
        }

        const typename Base::reference operator*() const override { return point_; }

        const Base& operator++() override {
            ++i_;
            if ( i_ == grid_.nx() ) {
                ++j_;
                i_ = 0;
            }
            compute_point_( i_, j_, point_ );
            return *this;
        }


        const Base& operator+=( typename Base::difference_type distance ) override {
            idx_t d = distance;
            while ( j_ != grid_.ny() && d >= ( grid_.nx() - i_ ) ) {
                d -= ( grid_.nx() - i_ );
                ++j_;
                i_ = 0;
            }
            i_ += d;
            compute_point_( i_, j_, point_ );
            return *this;
        }

        typename Base::difference_type distance( const Base& other ) const override {
            const auto& _other               = static_cast<const OrcaIterator&>( other );
            typename Base::difference_type d = 0;
            idx_t j                          = j_;
            idx_t i                          = i_;
            while ( j < _other.j_ ) {
                d += grid_.nx() - i;
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

    /// Constructor taking a path of an orca grid
    Orca( const std::string& name, const eckit::PathName& path );

private:
    /// Constructor taking a list of points (takes ownership)
    Orca( std::vector<PointXY>&& pts );

    /// Constructor taking a list of points (makes copy)
    Orca( const std::vector<PointXY>& pts );

    /// Constructor from initializer list
    Orca( std::initializer_list<PointXY> );

public:
    idx_t size() const override;

    Spec spec() const override;

    idx_t nx() const { return nx_; }
    idx_t ny() const { return ny_; }

    std::string name() const override;
    std::string type() const override;

    const PointXY& xy( idx_t i, idx_t j ) const { return points_[i + j * nx()]; }

    void xy( idx_t i, idx_t j, double crd[] ) const {
        const PointXY& p = xy( i, j );
        crd[0]           = p[0];
        crd[1]           = p[1];
    }

    PointLonLat lonlat( idx_t i, idx_t j ) const { return projection_.lonlat( xy( i, j ) ); }

    void lonlat( idx_t i, idx_t j, double crd[] ) const {
        xy( i, j, crd );
        projection_.xy2lonlat( crd );
    }

    virtual std::unique_ptr<Grid::IteratorXY> xy_begin() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this ) );
    }
    virtual std::unique_ptr<Grid::IteratorXY> xy_end() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this, false ) );
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_begin() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this ) );
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_end() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this, false ) );
    }

    size_t footprint() const override;

private:  // methods
    void print( std::ostream& ) const override;

    /// Hash of the lonlat array
    void hash( eckit::Hash& ) const override;

    RectangularLonLatDomain lonlatBoundingBox() const override;

private:
    /// Grid size
    idx_t nx_, ny_, periodicity_, halo_east_, halo_west_, halo_south_, halo_north_;

    /// Storage of coordinate points
    std::vector<PointXY> points_;

    /// name of the grid
    mutable std::string name_;

    /// Cache for the spec since may be quite heavy to compute
    mutable std::unique_ptr<Grid::Spec> cached_spec_;
};

extern "C" {
const Orca* atlas__grid__Orca__points( const double lonlat[], int shapef[], int stridesf[] );
const Orca* atlas__grid__Orca__config( util::Config* conf );
}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
