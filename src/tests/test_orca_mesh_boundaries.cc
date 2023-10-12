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
#include "atlas/output/Gmsh.h"

#include "atlas/util/Geometry.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {


CASE( "test haloExchange " ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        "eORCA1_T",  //
        "eORCA025_T",  //
    };
    for ( auto gridname : gridnames ) {
        for ( int64_t halo =0; halo < 2; ++halo ) {
            SECTION( gridname + "_halo" + std::to_string(halo) ) {
                auto grid = Grid(gridname);
                auto meshgen_config = grid.meshgenerator() | option::halo(halo);
                atlas::MeshGenerator meshgen(meshgen_config);
                auto partitioner_config = grid.partitioner();
                partitioner_config.set("type", "serial");
                auto partitioner = grid::Partitioner(partitioner_config);
                auto mesh = meshgen.generate(grid, partitioner);
                REQUIRE( mesh.grid() );
                EXPECT( mesh.grid().name() == gridname );
                idx_t count{0};

                const auto remote_idxs = array::make_indexview<idx_t, 1>(
                                          mesh.nodes().remote_index());
                functionspace::NodeColumns fs{mesh};
                Field field   = fs.createField<double>( option::name( "unswapped ghosts" ) );
                auto f        = array::make_view<double, 1>( field );
                const auto ghosts = atlas::array::make_view<int32_t, 1>(
                                      mesh.nodes().ghost());
                for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                    if (ghosts(jnode)) {
                        f( jnode ) = 1;
                    } else {
                        f( jnode ) = 0;
                    }
                }

                fs.haloExchange(field);

                for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                    if (f( jnode )) {
                      ++count;
                    }
                }
                if ( count != 0 ) {
                    Log::info() << "To diagnose problem, uncomment mesh writing here: " << Here() << std::endl;
                    //output::Gmsh gmsh(std::string("haloExchange_")+gridname+".msh",Config("coordinates","ij")|Config("info",true));
                    //gmsh.write(mesh);
                    //gmsh.write(field);
                }
                EXPECT_EQ( count, 0 );
            }
        }
    }
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
