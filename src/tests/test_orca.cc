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
    struct Section{
        std::string gridname;
        size_t size;
    };

    std::vector<Section> sections{
        {"ORCA2_T",27118},
        {"ORCA1_T",105704},
        {"ORCA025_T",1472282},
    };
    for( auto& section: sections ) {
        std::string gridname = section.gridname;
        SECTION( gridname ) {
            auto grid = Grid{gridname};

            EXPECT_EQ( grid.size(), section.size );

            Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

            idx_t n = 0;
            {
                auto trace = Trace(Here(),"iterating");
                for( auto& p : grid.lonlat() ) {
                    ++n;
                }
                trace.stop();
                Log::info() << "iterating took " << trace.elapsed() << " seconds" << std::endl;
            }
            EXPECT_EQ( n, grid.size() );
            Log::info() << "First point: " << grid.lonlat().front() << std::endl;
            Log::info() << "Last point: " << grid.lonlat().back() << std::endl;

            ATLAS_TRACE_SCOPE( "Mesh generation" ) {
                auto mesh = Mesh{grid};
                EXPECT_EQ(mesh.nodes().size(),grid.size());
            }

            // Now with extra virtual point at south pole
            EXPECT_EQ((MeshGenerator{"orca",util::Config("include_pole",true)}.generate(grid).nodes().size()),grid.size()+1);
            EXPECT_EQ((MeshGenerator{"orca",util::Config("force_include_south_pole",true)}.generate(grid).nodes().size()),grid.size()+1);
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
