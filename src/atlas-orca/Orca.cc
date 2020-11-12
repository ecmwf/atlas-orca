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

#include "eckit/log/ProgressTimer.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Statistics.h"
#include "eckit/utils/Hash.h"
#include "eckit/utils/Translator.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

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

    size_t npts = nx_halo_*ny_halo_;
    {
        eckit::Channel blackhole;
        eckit::ProgressTimer progress( "Reading "+name+" from file "+path_name.asString(), npts, " point", double( 1 ),
                                   npts > 2.e6 ? Log::info() : blackhole );

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
    if( trace.elapsed() > 1. ) {
        if( npts <= 2.e6 ){
            Log::info() << "Reading " << name << " from file " << path_name << " took " << trace.elapsed() << " seconds." << std::endl;
        }
        Log::info() << "  --> This will be greatly optimized soon when input is replaced with binary format." << std::endl;
    }
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

static class orca : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    orca() :
        GridBuilder( Orca::static_type(), {"^e?[Oo][Rr][Cc][Aa]([0-9]+)_([UVTF])$"}, {"[e]ORCA<deg>_{U,V,T,F}"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 30 ) << "[e]ORCA<deg>_{U,V,T,F}"
           << "ORCA Tripolar grid. Possible increasing resolutions <deg>: 2,1,025,12";
    }

    const Implementation* create( const std::string& name, const Config& /* config */ ) const override {
        int id;
        std::vector<std::string> matches;
        if ( not match( name, matches, id ) ) {
            return nullptr;
        }

        auto standard_name = []( std::string name ) {
            using atlas::orca::Library;
            auto to_upper = []( std::string str ) {
                std::for_each( str.begin(), str.end(), []( char& c ) {
                    c = static_cast<char>( std::toupper( static_cast<unsigned char>( c ) ) );
                } );
                return str;
            };
            name = to_upper( name );
            if ( name.front() == 'E' ) {
                name.front() = 'e';
            }
            return name;
        };
        auto computePath = []( std::string name ) {
            using atlas::orca::Library;
            auto to_lower = []( std::string str ) {
                std::for_each( str.begin(), str.end(), []( char& c ) {
                    c = static_cast<char>( std::tolower( static_cast<unsigned char>( c ) ) );
                } );
                return str;
            };
            name        = to_lower( name );
            name.back() = std::toupper( name.back() );
            return "~" + Library::instance().libraryName() + "/share/atlas-orca/data/" + name + ".ascii";
        };

        return new Orca{standard_name( name ), computePath( name )};
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
