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

#include "atlas-orca/Library.h"

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
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Spec.h"
#include "atlas/util/UnitSphere.h"


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

std::string spec_name( const Grid::Spec& spec, const std::string& def = "" ) {
    auto n = spec.getString( "orca_name", "" );
    auto a = spec.getString( "orca_arrangement", "" );
    return n.empty() || a.empty() ? def : ( n + "_" + a );
}

std::string spec_uid( const Grid::Spec& spec, const std::string& def = "" ) {
    return spec.getString( "uid", def );
}

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

Orca::Orca( const Config& config ) :
    Orca( spec_name( config, spec_uid( config, config.getString( "name", "" ) ) ), config ) {}

Orca::Orca( const std::string& name, const Config& config ) :
    name_( spec_name( config, spec_uid( config, name ) ) ), spec_( config ) {
#define EXPERIMENT_WITH_COARSENING 0

    auto trace = atlas::Trace( Here(), "Orca(" + name_ + ")" );

    eckit::PathName path = spec_.getString( "data" );
    if ( not path.exists() ) {
        ATLAS_THROW_EXCEPTION( "Could not locate orca grid file " << path );
    }

    std::ifstream ifstrm{path.asString().c_str()};

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
        eckit::ProgressTimer progress( "Reading " + name_ + " from file " + path, npts, " point", double( 1 ),
                                       npts > 2.e6 ? Log::info() : blackhole );

        PointXY p;
        double lsm, core;
        for ( auto jj = 0; jj != ny_halo_; ++jj ) {
            for ( auto ii = 0; ii != nx_halo_; ++ii, ++progress ) {
                std::getline( ifstrm, line );
                std::istringstream iss{line};
                ATLAS_ASSERT( iss >> p[1] >> p[0] >> lsm >> core, "Error while reading coordinates" );
                lsm_.push_back( lsm );
                core_.push_back( core );
                points_halo_.emplace_back( p );
            }
        }
    }

    trace.stop();
    if ( trace.elapsed() > 1. ) {
        if ( npts <= 2.e6 ) {
            Log::info() << "Reading " << name_ << " from file " << path << " took " << trace.elapsed() << " seconds."
                        << std::endl;
        }
        Log::info() << "  --> This will be greatly optimized soon when input is replaced with binary format."
                    << std::endl;
    }


    // Fixes, TBD
    // If these fixes are not applied, halo exchanges or other numerics may not be performed correctly
    double west = 0.;

    {
        auto point    = [&]( idx_t i, idx_t j ) -> PointXY& { return const_cast<PointXY&>( xy( i, j ) ); };
        auto set_core = [&]( idx_t i, idx_t j, bool core ) { core_[( imin_ + i ) + ( jmin_ + j ) * jstride_] = core; };

        if ( name_ == "ORCA12_T" ) {
            set_core( 2160, 3056, true );
        }
        if ( name_ == "eORCA12_T" ) {
            set_core( 2160, 3603, true );
        }
        if ( name_ == "ORCA025_T" ) {
            set_core( 720, 1018, true );
            point( 720, 1019 ) = xy( 720, 1017 );
        }
        if ( name_ == "eORCA025_T" ) {
            set_core( 720, 1204, true );
        }
        if ( name_ == "ORCA2_T" ) {
            set_core( 90, 146, true );
        }
        if ( name_ == "ORCA1_T" ) {
            point( -1, 290 )  = xy( 0, 289 );
            point( 359, 290 ) = xy( 0, 289 );
            point( 360, 290 ) = xy( 359, 289 );
        }
        if ( name_ == "eORCA1_T" ) {
            point( -1, 330 )  = xy( 0, 329 );
            point( 359, 330 ) = xy( 0, 329 );
            point( 360, 330 ) = xy( 359, 329 );
        }
        if ( name_ == "ORCA1_F" ) {
            point( -1, 290 )  = xy( 359, 289 );
            point( 359, 290 ) = xy( 359, 289 );
        }

        west = spec_.getString( "orca_name" ) == "ORCA2" ? 80. : 73.;
    }

    domain_ = GlobalDomain( west );

    info_.reset( new OrcaInfo( *this ) );
}

void Orca::hash( eckit::Hash& h ) const {
    h.add( spec_uid( spec_ ) );
}

idx_t Orca::size() const {
    return nx_halo_ * ny_halo_;
}

/// @return parallel/meridian limits containing the grid
RectangularLonLatDomain Orca::lonlatBoundingBox() const {
    return domain_;
}

Grid::Spec Orca::spec() const {
    return spec_;
}

void Orca::print( std::ostream& os ) const {
    os << "Orca(" << name_ << ")";
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

static class OrcaGridBuilder : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    OrcaGridBuilder() : GridBuilder( Orca::static_type(), {"^e?ORCA[0-9]+_[FTUVW]$"}, {"[e]ORCA<deg>_{F,T,U,V,W}"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 30 ) << "[e]ORCA<deg>_{F,T,U,V,W}"
           << "ORCA Tripolar grid. Possible increasing resolutions <deg>: 2,1,025,12";
    }

    const Implementation* create( const std::string& name_or_uid, const Config& /* config */ ) const override {
        using Registry = util::SpecRegistry<atlas::Grid>;

        auto sane_id( name_or_uid );
        std::transform( sane_id.begin(), sane_id.end(), sane_id.begin(), ::tolower );

        if ( Registry::has( sane_id ) ) {
            return create( Registry::lookup( sane_id ) );
        }

        auto sane_name( name_or_uid );
        std::transform( sane_name.begin(), sane_name.end(), sane_name.begin(), ::toupper );
        if ( sane_name.front() == 'E' ) {
            sane_name.front() = 'e';
        }

        if ( Registry::has( sane_name ) ) {
            return create( Registry::lookup( sane_name ) );
        }

        return nullptr;
    }

    const Implementation* create( const Config& config ) const override { return new Orca( config ); }

    void force_link() {}

} orca_;

void force_link_Orca() {
    orca_.force_link();
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
