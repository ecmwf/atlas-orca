/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Orca.h"
#include "Library.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>

#include "atlas/domain/Domain.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"
#include "atlas/util/UnitSphere.h"

#include "eckit/filesystem/PathName.h"
#include "eckit/log/Statistics.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

std::string Orca::static_type() {
    return "orca";
}

std::string Orca::name() const {
    return name_;
}

std::string Orca::type() const {
    return static_type();
}

Orca::Orca( const std::string& name, const eckit::PathName& path_name ) : name_( name ) {
#define EXPERIMENT_WITH_COARSENING 0

    ATLAS_TRACE( "Orca(" + name + ")" );

    if ( not path_name.exists() ) {
        ATLAS_THROW_EXCEPTION( "Could not locate orca grid file " << path_name );
    }

    std::ifstream ifstrm{path_name.asString().c_str()};

    // Reading header
    std::string line;
    std::getline( ifstrm, line );
    std::istringstream iss{line};

    idx_t nx, ny;
    ATLAS_ASSERT( iss >> nx >> ny >> periodicity_ >> halo_east_ >> halo_west_ >> halo_south_ >> halo_north_,
                  "Error while reading header " );

    nx_ = nx - halo_east_ - halo_west_;
    ny_ = ny - halo_north_ - halo_south_;
    points_.reserve( nx_ * ny_ );

#if EXPERIMENT_WITH_COARSENING
    // halo_south_ = 9 * ny / 10;   // Increase to decrease latitudes from South Pole
    int skipx = 2;
    int skipy = 2;
    nx_       = 0;
    ny_       = 0;
#endif
    auto is_halo = [&]( idx_t i, idx_t j ) {
        if ( i < halo_west_ ) {
            return true;
        }
        if ( i >= nx - halo_east_ ) {
            return true;
        }
        if ( j < halo_south_ ) {
            return true;
        }
        if ( j >= ny - halo_north_ ) {
            return true;
        }
        return false;
    };

    // Reading coordinates
    PointXY p;
    for ( auto jj = 0; jj != ny; ++jj ) {
#if EXPERIMENT_WITH_COARSENING
        bool consider_j     = false;
        idx_t diff_from_top = ny - halo_north_ - 1 - jj;
        if ( diff_from_top % skipy == 0 ) {
            consider_j = true;
        }
        bool j_added = false;
        idx_t _nx    = 0;
#endif

        for ( auto ii = 0; ii != nx; ++ii ) {
            std::getline( ifstrm, line );
            std::istringstream iss{line};
            ATLAS_ASSERT( iss >> p[1] >> p[0], "Error while reading coordinates" );

            if ( is_halo( ii, jj ) ) {
                continue;
            }
#if !EXPERIMENT_WITH_COARSENING
            points_.emplace_back( p );
        }
#else
            if ( consider_j ) {
                j_added = true;
                if ( ( halo_west_ + ii ) % skipx == 0 ) {
                    points_.push_back( p );
                    ++_nx;
                }
            }
        }
        if ( jj == ny - halo_north_ - 1 ) {
            ATLAS_ASSERT( j_added );
        }
        if ( j_added ) {
            if ( nx_ == 0 ) {
                nx_ = _nx;
            }
            ++ny_;
        }
#endif
    }
    ATLAS_ASSERT( points_.size() == nx_ * ny_ );

    domain_ = GlobalDomain{};
}

void Orca::hash( eckit::Hash& h ) const {}

idx_t Orca::size() const {
    return static_cast<idx_t>( points_.size() );
}

/// @return parallel/meridian limits containing the grid
RectangularLonLatDomain Orca::lonlatBoundingBox() const {
    return projection_ ? projection_.lonlatBoundingBox( domain_ ) : domain_;
}

Grid::Spec Orca::spec() const {
    return util::Config{"name", name()};
}

void Orca::print( std::ostream& os ) const {
    os << "Orca(" << name() << ")";
}

size_t Orca::footprint() const {
    return sizeof( PointXY ) * points_.size();
}

namespace {

static class orca : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    orca() : GridBuilder( Orca::static_type(), {"^e?orca([0-9]+)_([UVTF])$"}, {"[e]orca<deg>_{U,V,T,F}"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 30 ) << "[e]orca<deg>_{U,V,T,F}"
           << "ORCA Tripolar grid. Possible increasing resolutions <deg>: 2,1,025,12";
    }

    const Implementation* create( const std::string& name, const Config& /* config */ ) const override {
        int id;
        std::vector<std::string> matches;
        if ( not match( name, matches, id ) ) {
            return nullptr;
        }

        auto computePath = []( std::string name ) {
            using atlas::orca::Library;
            return "~" + Library::instance().libraryName() + "/share/atlas-orca/orca/" + name + ".ascii";
        };

        return new Orca{name, computePath( name )};
    }

    const Implementation* create( const Config& config ) const override {
        throw_NotImplemented( "Cannot create unstructured grid from config", Here() );
    }

    void force_link() {}

} orca_;
}  // namespace

void force_link_Orca() {
    orca_.force_link();
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
