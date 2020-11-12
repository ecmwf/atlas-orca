/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/system/ResourceUsage.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "atlas-orca/OrcaGrid.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

size_t peakMemory() {
    return eckit::system::ResourceUsage().maxResidentSetSize();
}

CASE( "test generate orca mesh" ) {
    std::string gridname = "ORCA1_V";
    SECTION( gridname ) {
        auto grid = Grid{gridname};
        Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

        auto meshgenerator = MeshGenerator{"orca"};
        auto mesh          = meshgenerator.generate( grid );
        Log::info() << "mesh.footprint() = " << eckit::Bytes( mesh.footprint() ) << std::endl;

        if ( mesh.footprint() < 25 * 1e6 ) {  // less than 25 Mb
            output::Gmsh{"orca_2d.msh", Config( "coordinates", "lonlat" )}.write( mesh );
            output::Gmsh{"orca_3d.msh", Config( "coordinates", "xyz" )}.write( mesh );
        }
        ATLAS_DEBUG( "Peak memory: " << eckit::Bytes( peakMemory() ) );
    }

    SECTION( "auto_generate" ) { auto mesh = Mesh{Grid{gridname}}; }
}


CASE( "test orca grid iterator" ) {
    std::string gridname = "ORCA2_T";
    SECTION( gridname ) {
        auto grid = Grid{gridname};
        Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

        idx_t n = 0;
        ATLAS_TRACE_SCOPE("iterating") {
            for( auto& p : grid.lonlat() ) {
                ++n;
            }
        }
        EXPECT_EQ( n, grid.size() );
        Log::info() << "First point: " << grid.lonlat().front() << std::endl;
        Log::info() << "Last point: " << grid.lonlat().back() << std::endl;
        auto mesh = Mesh{grid};
        EXPECT_EQ(mesh.nodes().size(),grid.size()+1); // One extra virtual point at south pole
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
