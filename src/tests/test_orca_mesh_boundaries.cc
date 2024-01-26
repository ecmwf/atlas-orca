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
#include <sstream>

#include "atlas-orca/meshgenerator/OrcaMeshGenerator.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/util/function/VortexRollup.h"

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


CASE( "test haloExchange " ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        "ORCA2_U",   //
        "ORCA2_V",   //
        "eORCA1_T",  //
        "eORCA025_T",  //
        "eORCA12_T",  //
    };
    auto distributionNames = std::vector<std::string>{
      "serial",
      "checkerboard",
      "equal_regions",
      // "equal_area",
    };

    auto rollup_plus = []( const double lon, const double lat ) {
        return 1 + util::function::vortex_rollup( lon, lat, 0.0 );
    };

    for ( auto distributionName : distributionNames ) {
        for ( auto gridname : gridnames ) {
            for ( int64_t halo = 0; halo < 1; ++halo ) {
                if ( distributionName == "serial" && halo != 0 )
                  continue;
                SECTION( gridname + "_" + distributionName + "_halo" + std::to_string(halo) ) {
                    auto grid = Grid(gridname);
                    auto meshgen_config = grid.meshgenerator() | option::halo(halo);
                    atlas::MeshGenerator meshgen(meshgen_config);
                    auto partitioner_config = grid.partitioner();
                    partitioner_config.set( "type", distributionName );
                    auto partitioner = grid::Partitioner( partitioner_config );
                    auto mesh        = meshgen.generate( grid, partitioner );
                    std::cout << "mesh generator finished " << std::endl;
                    REQUIRE( mesh.grid() );
                    EXPECT( mesh.grid().name() == gridname );
                    idx_t count{ 0 };

                    const auto remote_idxs = array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() );
                    functionspace::NodeColumns fs( mesh, option::halo( halo ) );
                    double halosize = 0;
                    mesh.metadata().get( "halo", halosize );
                    EXPECT( halosize == halo );

                    Field field       = fs.createField<double>( option::name( "vortex rollup" ) );
                    Field field2      = fs.createField<double>( option::name( "remotes < 0" ) );
                    auto f            = array::make_view<double, 1>( field );
                    auto f2           = array::make_view<double, 1>( field2 );
                    const auto ghosts = atlas::array::make_view<int32_t, 1>( mesh.nodes().ghost() );
                    const auto lonlat = array::make_view<double, 2>( mesh.nodes().lonlat() );
                    std::cout << "begin loop over mesh nodes " << std::endl;
                    for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                        if ( ghosts( jnode ) ) {
                            f( jnode ) = 0;
                        }
                        else {
                            const double lon = lonlat( jnode, 0 );
                            const double lat = lonlat( jnode, 1 );
                            f( jnode )       = rollup_plus( lon, lat );
                        }
                    }

                    fs.haloExchange( field );

                    const auto xy        = array::make_view<double, 2>( mesh.nodes().xy() );
                    const auto glb_idxs  = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );
                    const auto partition = array::make_view<int32_t, 1>( mesh.nodes().partition() );
                    const auto halos     = array::make_view<int32_t, 1>( mesh.nodes().halo() );

                    const auto master_glb_idxs =
                        array::make_view<gidx_t, 1>( mesh.nodes().field( "master_global_index" ) );
                    const auto ij = array::make_view<idx_t, 2>( mesh.nodes().field( "ij" ) );

                    double sumSquares{ 0.0 };
                    int halocount = 0;
                    for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                        f2( jnode )      = 0;
                        const double lon = lonlat( jnode, 0 );
                        const double lat = lonlat( jnode, 1 );
                        if ( std::abs( f( jnode ) - rollup_plus( lon, lat ) ) > 1e-6 ) {
                            f( jnode ) = -1;
                            sumSquares += std::pow( std::abs( f( jnode ) - rollup_plus( lon, lat ) ), 2 );
                            ++count;
                        }
                        if ( remote_idxs( jnode ) < 0 ) {
                            std::cout << "[" << mpi::rank() << "] remote_idx < 0 " << jnode << " : ghost "
                                      << ghosts( jnode ) << " master_global_index " << master_glb_idxs( jnode )
                                      << " lon " << lonlat( jnode, 0 ) << " lat " << lonlat( jnode, 1 ) << std::endl;
                            f2( jnode ) = 1;
                        }
                        if (halos(jnode) > 0) {
                          std::cout << "[" << mpi::rank() << "] i " << ij(jnode, 0) << " j " << ij(jnode, 1)
                                    << " halo(" << jnode << ") " << halos(jnode)
                                    << " ghost " << ghosts(jnode) << " master_global_index " << master_glb_idxs(jnode)
                                    << " lon " << lonlat(jnode, 0) << " lat " << lonlat(jnode, 1) << std::endl;
                          halocount++;
                        }
                    }
                    //if ( count != 0 ) {
                        Log::info() << "count nonzero and norm of differences is: " << std::sqrt(sumSquares) << std::endl;
                        Log::info() << "count nonzero halo points: " << halocount << std::endl;
                        Log::info() << "To diagnose problem, uncomment mesh writing here: " << Here() << std::endl;
                        //output::Gmsh gmsh(
                        //    std::string("haloExchange_")+gridname+"_"+distributionName+"_"+std::to_string(halo)+".msh",
                        //    Config("coordinates","ij")|Config("info",true));
                        //gmsh.write(mesh);
                        //gmsh.write(field);
                        //gmsh.write(field2);
                    //}
                    EXPECT_EQ( count, 0 );
                }
            }
        }
    }
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
