/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas-orca/Orca.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>

#include "eckit/filesystem/PathName.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/log/Statistics.h"
#include "eckit/utils/Hash.h"
#include "eckit/utils/Translator.h"

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
#include "atlas/util/Spec.h"
#include "atlas/util/UnitSphere.h"

#include "atlas-orca/Library.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {


namespace {

struct PointIJ {
    idx_t i;
    idx_t j;
    friend std::ostream& operator<<( std::ostream& out, const PointIJ& p ) {
        out << "{" << p.j << "," << p.i << "}";
        return out;
    }
    bool operator==( const PointIJ& other ) const { return other.i == i && other.j == j; }
    bool operator!=( const PointIJ& other ) const { return other.i != i || other.j != j; }
    bool operator<( const PointIJ& other ) const {
        if ( j < other.j ) {
            return true;
        }
        if ( j == other.j && i < other.i ) {
            return true;
        }
        return false;
    }
};
}  // namespace
struct Orca::OrcaInfo {
    OrcaInfo( const Orca& g ) : grid( g ) {
        ix_pivot = grid.nx() / 2;
        if ( g.name() == "ORCA1_F" ) {
            ix_pivot--;
            exceptions[PointIJ{-1, 290}]  = PointIJ{359, 289};
            exceptions[PointIJ{359, 290}] = PointIJ{359, 289};
        }
        patch = not grid.ghost( ix_pivot + 1, grid.ny() - 1 );
        //        ATLAS_DEBUG_VAR( g.name() );
        //        ATLAS_DEBUG_VAR( g.nx() );
        //        ATLAS_DEBUG_VAR( g.ny() );
        //        ATLAS_DEBUG_VAR( g.haloWest() );
        //        ATLAS_DEBUG_VAR( g.haloNorth() );
        //        ATLAS_DEBUG_VAR( patch );
        //        ATLAS_DEBUG_VAR( ix_pivot );
    }
    bool patch;
    int ix_pivot;
    std::map<PointIJ, PointIJ> exceptions;
    idx_t i_valid( idx_t i ) const {
        if ( i < 0 ) {
            return i + grid.nx();
        }
        else if ( i >= grid.nx() ) {
            return i - grid.nx();
        }
        else {
            return i;
        }
    }
    bool in_halo( idx_t i ) const { return i < 0 || i >= grid.nx(); }
    PointIJ periodicIndex( idx_t i, idx_t j ) const {
        if ( exceptions.size() ) {
            auto it = exceptions.find( PointIJ{i, j} );
            if ( it != exceptions.end() ) {
                return it->second;
            }
        }
        PointIJ p;
        if ( not patch ) {  // e.g. ORCA2_T, ORCA025_T, ORCA12_T, ORCA1_F
            int iy_pivot = grid.ny() - 1;
            if ( grid.ghost( i, j ) ) {
                i = i_valid( i );
                if ( grid.ghost( i, j ) && j >= 0 ) {
                    j = iy_pivot - ( j - iy_pivot );
                    i = ix_pivot - ( i - ix_pivot );
                    i = i_valid( i );
                }
                p.i = i;
                p.j = j;
            }
            else {
                p.j = iy_pivot - ( j - iy_pivot );
                if ( i < 0 || ( i <= 0 && j >= grid.ny() ) ) {
                    p.i = -i;
                }
                else {
                    p.i = ix_pivot - ( i - ix_pivot );
                }
            }
        }
        else {  // patch, e.g. ORCA1_T
            Log::Channel hole;
            double iy_pivot = double( grid.ny() ) - 0.5;

            if ( grid.ghost( i, j ) ) {
                i = i_valid( i );
                if ( grid.ghost( i, j ) && j >= 0 ) {
                    i = 2 * ix_pivot - 1 - i;
                    j = int( 2. * iy_pivot ) - j;
                    i = i_valid( i );
                }
                p.i = i;
                p.j = j;
            }
            else {
                p.j = int( 2. * iy_pivot ) - j;

                if ( i < 0 && j >= grid.ny() ) {
                    p.i = -i - 1;
                }
                else if ( i <= 0 && j < grid.ny() ) {
                    p.i = i + grid.nx();
                    p.j = j;
                }
                else if ( i >= grid.nx() && j < grid.ny() ) {
                    p.i = i - grid.nx();
                    p.j = j;
                }
                else if ( i >= grid.nx() && j >= grid.ny() ) {
                    p.i = i - grid.nx();
                }
                else {
                    p.i = 2 * ix_pivot - 1 - i;
                }
            }
        }
        if ( p.i < -grid.haloWest() ) {
            p.i += grid.nx();
        }
        return p;
    }

private:
    const Orca& grid;
};

namespace {
double hardcoded_west( const std::string& name ) {  // should come from file instead!
    static std::map<std::string, double> west = []() {
        std::map<std::string, double> west;
        std::vector<std::string> prefixes = {"ORCA", "eORCA"};
        std::vector<std::string> suffixes = {"_T", "_U", "_V", "_F"};
        for ( auto& prefix : prefixes ) {
            for ( auto& suffix : suffixes ) {
                west[prefix + "2" + suffix]   = 80.;
                west[prefix + "1" + suffix]   = 73.;
                west[prefix + "025" + suffix] = 73.;
                west[prefix + "12" + suffix]  = 73.;
            }
        }
        return west;
    }();
    return west.at( name );
}

}  // namespace


std::string Orca::static_type() {
    return "ORCA";
}

std::string Orca::name() const {
    return name_;
}

std::string Orca::type() const {
    return static_type();
}

gidx_t Orca::periodicIndex( idx_t i, idx_t j ) const {
    // Something like this will not need to be computed in the future, it will come from file.
    PointIJ p = info_->periodicIndex( i, j );
    return index( p.i, p.j );
}

Orca::Orca( const std::string& name, const eckit::PathName& path_name ) : name_( name ) {
#define EXPERIMENT_WITH_COARSENING 0

    auto trace = atlas::Trace( Here(), "Orca(" + name + ")" );

    if ( not path_name.exists() ) {
        ATLAS_THROW_EXCEPTION( "Could not locate orca grid file " << path_name );
    }

    std::ifstream ifstrm{path_name.asString().c_str()};

    // Reading header
    std::string line;
    std::getline( ifstrm, line );
    std::istringstream iss{line};

    ATLAS_ASSERT( iss >> nx_halo_ >> ny_halo_ >> periodicity_ >> halo_east_ >> halo_west_ >> halo_south_ >> halo_north_,
                  "Error while reading header " );

    auto getEnv = []( const std::string& env, int default_value ) {
        if ( ::getenv( env.c_str() ) ) {
            return eckit::Translator<std::string, int>()( ::getenv( env.c_str() ) );
        }
        return default_value;
    };

    halo_south_ = getEnv( "halo_south", halo_south_ );

    nx_      = nx_halo_ - halo_east_ - halo_west_;
    ny_      = ny_halo_ - halo_north_ - halo_south_;
    jstride_ = nx_halo_;

    halo_north_ = std::min( halo_north_, getEnv( "halo_north", halo_north_ ) );
    imin_       = halo_east_;
    jmin_       = halo_south_;

    points_halo_.reserve( nx_halo_ * ny_halo_ );
    lsm_.reserve( nx_halo_ * ny_halo_ );
    core_.reserve( nx_halo_ * ny_halo_ );

    // Reading coordinates

    size_t npts = nx_halo_ * ny_halo_;
    {
        eckit::Channel blackhole;
        eckit::ProgressTimer progress( "Reading " + name + " from file " + path_name.asString(), npts, " point",
                                       double( 1 ), npts > 2.e6 ? Log::info() : blackhole );

        PointXY p;
        double lsm, core;
        for ( auto jj = 0; jj != ny_halo_; ++jj ) {
            for ( auto ii = 0; ii != nx_halo_; ++ii, ++progress ) {
                std::getline( ifstrm, line );
                std::istringstream iss{line};
                ATLAS_ASSERT( iss >> p[1] >> p[0] >> lsm >> core, "Error while reading coordinates" );
                lsm_.emplace_back( lsm );
                core_.emplace_back( core );
                points_halo_.emplace_back( p );
            }
        }
    }
    domain_ = GlobalDomain{hardcoded_west( name )};
    trace.stop();
    if ( trace.elapsed() > 1. ) {
        if ( npts <= 2.e6 ) {
            Log::info() << "Reading " << name << " from file " << path_name << " took " << trace.elapsed()
                        << " seconds." << std::endl;
        }
        Log::info() << "  --> This will be greatly optimized soon when input is replaced with binary format."
                    << std::endl;
    }


    // Fixes, TBD
    // If these fixes are not applied, halo exchanges or other numerics may not be performed correctly
    {
        auto point    = [&]( idx_t i, idx_t j ) -> PointXY& { return const_cast<PointXY&>( xy( i, j ) ); };
        auto set_core = [&]( idx_t i, idx_t j, bool core ) { core_[( imin_ + i ) + ( jmin_ + j ) * jstride_] = core; };
        if ( name == "ORCA12_T" ) {
            set_core( 2160, 3056, true );
        }
        if ( name == "ORCA025_T" ) {
            set_core( 720, 1018, true );
            point( 720, 1019 ) = xy( 720, 1017 );
        }
        if ( name == "ORCA2_T" ) {
            set_core( 90, 146, true );
        }
        if ( name == "ORCA1_T" ) {
            point( -1, 290 )  = xy( 0, 289 );
            point( 359, 290 ) = xy( 0, 289 );
            point( 360, 290 ) = xy( 359, 289 );
        }
        if ( name == "ORCA1_F" ) {
            point( -1, 290 )  = xy( 359, 289 );
            point( 359, 290 ) = xy( 359, 289 );
        }
    }

    info_.reset( new OrcaInfo( *this ) );
}

void Orca::hash( eckit::Hash& h ) const {
    h.add( name_ );
}

idx_t Orca::size() const {
    return nx_halo_ * ny_halo_;
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

namespace {
template <typename T>
size_t memory( const T& container ) {
    return sizeof( typename T::value_type ) * container.size();
}
}  // namespace
size_t Orca::footprint() const {
    return memory( points_halo_ ) + memory( lsm_ ) + memory( core_ );
}

Grid::Config Orca::meshgenerator() const {
    return Config( "type", "orca" );
}

Grid::Config Orca::partitioner() const {
    return Config( "type", "checkerboard" );
}

namespace {

static class OrcaGridBuilder : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    OrcaGridBuilder() :
        GridBuilder( Orca::static_type(), {"^e?[Oo][Rr][Cc][Aa]([0-9]+)_([UVTF])$"}, {"[e]ORCA<deg>_{U,V,T,F}"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 30 ) << "[e]ORCA<deg>_{U,V,T,F}"
           << "ORCA Tripolar grid. Possible increasing resolutions <deg>: 2,1,025,12";
    }

    const Implementation* create( const std::string& name, const Config& /* config */ ) const override {
        {
            int id;
            std::vector<std::string> matches;
            if ( not match( name, matches, id ) ) {
                return nullptr;
            }
        }

        auto id = name;
        std::transform( id.begin(), id.end(), id.begin(), ::toupper );
        if ( id.front() == 'E' ) {
            id.front() = 'e';
        }

        return create( util::SpecRegistry<atlas::Grid>::lookup( id ) );
    }

    const Implementation* create( const Config& config ) const override {
        return new Orca( config.getString( "name", "" ), config.getString( "data" ) );
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
