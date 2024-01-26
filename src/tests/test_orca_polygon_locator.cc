/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <numeric>

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Config.h"

#include "atlas/util/PolygonLocator.h"
#include "atlas/util/PolygonXY.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test orca polygon locator" ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",     //
        "eORCA1_T",    //
        "eORCA025_T",  //
        "eORCA12_T",  //
    };

    std::string grid_resource = eckit::Resource<std::string>( "--grid", "" );
    if ( not grid_resource.empty() ) {
        gridnames = { grid_resource };
    }

    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            OrcaGrid grid = Grid( gridname );
            grid::Partitioner partitioner( "equal_regions", atlas::mpi::size() );
            auto meshgenerator = MeshGenerator{ "orca" };
            auto mesh          = meshgenerator.generate( grid, partitioner );
            REQUIRE( mesh.grid() );
            EXPECT( mesh.grid().name() == gridname );
            atlas::util::ListPolygonXY polygon_list( mesh.polygons() );
            //mesh.polygon(0).outputPythonScript("polygon.py",
            //                                   Config("nodes", false)("coordinates", "xy"));
            atlas::util::PolygonLocator locator( polygon_list, mesh.projection() );

            // fails on ORCA2 because partition polygon is not connected space
            // (ORCA2 grids have a cut out area for the mediterranean)
            idx_t part;
            part = locator( { 82.0, 10.0 } );
            ATLAS_DEBUG_VAR( part );
            // would fail on ORCA025 because xy coordinate system doesn't cover
            // 0-360. A fix for this has been added to the locator in atlas 0.32.0.
            part = locator( { -173.767, -61.1718 } );
            ATLAS_DEBUG_VAR( part );
        }
    }
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
