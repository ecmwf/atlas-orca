/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/system/ResourceUsage.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "atlas/util/Geometry.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

size_t peakMemory() {
    return eckit::system::ResourceUsage().maxResidentSetSize();
}

//-----------------------------------------------------------------------------

CASE( "test generate orca mesh" ) {
    static std::string gridname = "eORCA1_T";
    static auto grid            = Grid( gridname );

    SECTION( "orca_generate" ) {
        Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

        auto meshgenerator = MeshGenerator{"orca"};
        auto mesh          = meshgenerator.generate( grid );
        Log::info() << "mesh.footprint() = " << eckit::Bytes( mesh.footprint() ) << std::endl;

        EXPECT_EQ( mesh.nodes().size(), grid.size() );

        if ( mesh.footprint() < 25 * 1e6 ) {  // less than 25 Mb
            output::Gmsh{"orca_2d.msh", Config( "coordinates", "lonlat" )}.write( mesh );
            output::Gmsh{"orca_3d.msh", Config( "coordinates", "xyz" )}.write( mesh );
        }
        ATLAS_DEBUG( "Peak memory: " << eckit::Bytes( peakMemory() ) );
    }

    SECTION( "auto_generate" ) { auto mesh = Mesh{grid}; }
}

//-----------------------------------------------------------------------------

CASE( "test orca mesh halo" ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        "eORCA1_T",  //
        //"eORCA025_T",  //
    };
    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            auto mesh = Mesh{Grid( gridname )};
            REQUIRE( mesh.grid() );
            EXPECT( mesh.grid().name() == gridname );
            auto remote_idx = array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() );
            auto ij         = array::make_view<idx_t, 2>( mesh.nodes().field( "ij" ) );
            idx_t count{0};

            functionspace::NodeColumns fs{mesh};
            Field field   = fs.createField<double>( option::name( "bla" ) );
            auto f        = array::make_view<double, 1>( field );
            OrcaGrid grid = mesh.grid();
            for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                if ( remote_idx( jnode ) < 0 ) {
                    auto p = orca::PointIJ{ij( jnode, 0 ), ij( jnode, 1 )};
                    orca::PointIJ master;
                    grid->index2ij( grid->periodicIndex( p.i, p.j ), master.i, master.j );
                    Log::info() << p << " --> " << master << std::endl;
                    ++count;
                    f( jnode ) = 1;
                }
                else {
                    f( jnode ) = 0;
                }
            }
            EXPECT_EQ( count, 0 );
            if ( count != 0 ) {
                Log::info() << "To diagnose problem, uncomment mesh writing here: " << Here() << std::endl;
                //                output::Gmsh gmsh(gridname+".msh",Config("coordinates","ij")|Config("ghost",true)|Config("info",true));
                //                gmsh.write(mesh);
                //                gmsh.write(field);
            }
        }
    }
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
