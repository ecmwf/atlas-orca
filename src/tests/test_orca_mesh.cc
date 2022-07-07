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
#include "atlas/util/Unique.h"

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

void build_remote_halo_info(Mesh& mesh, array::ArrayView<int, 1> remote_halo );
void build_periodic_boundaries( Mesh& mesh );  // definition below

//-----------------------------------------------------------------------------

CASE( "test generate orca mesh" ) {
    static std::string gridname = "ORCA2_T";
    static auto grid            = Grid( gridname );

    SECTION( "orca_generate" ) {
        Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

        auto meshgenerator = MeshGenerator{"orca", option::halo( 1 ) };
        auto mesh          = meshgenerator.generate( grid );
        Log::info() << "mesh.footprint() = " << eckit::Bytes( mesh.footprint() ) << std::endl;

        if (mpi::comm().size() == 1) EXPECT_EQ( mesh.nodes().size(), grid.size() );

        if ( mesh.footprint() < 25 * 1e6 ) {  // less than 25 Mb
            output::Gmsh{"2d.msh", Config( "coordinates", "lonlat" )}.write( mesh );
            output::Gmsh{"3d.msh", Config( "coordinates", "xyz" )}.write( mesh );
        }
        ATLAS_DEBUG( "Peak memory: " << eckit::Bytes( peakMemory() ) );
    }

    SECTION( "auto_generate" ) { auto mesh = Mesh{grid}; }
}

//-----------------------------------------------------------------------------

CASE( "test orca mesh halo" ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        //"eORCA1_T",  //
        //"eORCA025_T",  //
        //"eORCA12_T",  //
    };
    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            OrcaGrid grid      = Grid( gridname );
            grid::Partitioner partitioner("equal_regions", atlas::mpi::size());
            auto meshgenerator = MeshGenerator{"orca", option::halo( 2 ) };
            auto mesh          = meshgenerator.generate( grid, partitioner );
            REQUIRE( mesh.grid() );
            EXPECT( mesh.grid().name() == gridname );
            int mpi_rank = mpi::rank();
            auto remote_idx = array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() );
            auto ij         = array::make_view<idx_t, 2>( mesh.nodes().field( "ij" ) );
            idx_t count{0};

            functionspace::NodeColumns fs{mesh,  option::halo( 2 )};
            Field field   = fs.createField<int>( option::name( "remote_halo" ));
            auto f        = array::make_view<int, 1>( field );
            auto master_glb_index = array::make_view<gidx_t, 1>(mesh.nodes().field("master_global_index"));
            auto glb_index = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
            for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                f(jnode) = 0;
                if ( remote_idx( jnode ) < 0 ) {
                    auto p = orca::PointIJ{ij( jnode, 0 ), ij( jnode, 1 )};
                    orca::PointIJ master;
                    grid->index2ij( grid->periodicIndex( p.i, p.j ), master.i, master.j );
                    //std::cout << "[" << mpi_rank << "] " << p << " --> " << master << std::endl;
                    ++count;
                }
            }
            fs.halo_exchange();
            if ( count != 0 ) {
                build_remote_halo_info( mesh, f );
                Log::info() << "To diagnose problem, uncomment mesh writing here: " << Here() << std::endl;
                output::Gmsh gmsh(gridname+".msh",Config("coordinates","ij")|Config("ghost",true)|Config("info",true));
                gmsh.write(mesh);
                gmsh.write(field);
            }
            EXPECT_EQ( count, 0 );
        ATLAS_DEBUG( "Peak memory: " << eckit::Bytes( peakMemory() ) );
        }
    }
}

using Unique2Node = std::map<gidx_t, idx_t>;
void build_remote_halo_info(Mesh& mesh, array::ArrayView<int, 1> remote_halo ) {
    ATLAS_TRACE();
    auto mpi_size = mpi::size();
    auto mypart   = mpi::rank();
    mesh::Nodes& nodes = mesh.nodes();
    int nb_nodes = nodes.size();

    // get hold of the halo data remote indices, partitions, and halo values
    auto gidx    = array::make_indexview<gidx_t, 1>(nodes.field("master_global_index"));
    auto ridx    = array::make_indexview<idx_t, 1>( nodes.remote_index() );
    auto part    = array::make_view<int, 1>( nodes.partition() );
    auto ghost   = array::make_view<int, 1>( nodes.ghost() );
    auto halo    = array::make_view<int, 1>( nodes.halo() );

    // organise the halo data into send buffers
    std::vector<std::vector<int>> send_halo( mpi_size );
    std::vector<std::vector<int>> send_gidx( mpi_size );

    Unique2Node global2local;

    for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
        gidx_t uid     = gidx(jnode);
        if (halo (jnode) != 0) {
            send_halo[part(jnode)].push_back(halo(jnode));
            send_gidx[part(jnode)].push_back(uid);
        }
        if (not ghost(jnode)) {
            bool inserted = global2local.insert(std::make_pair(uid, jnode)).second;
            ATLAS_ASSERT(inserted, std::string("index already inserted ") + std::to_string(uid) + ", "
                + std::to_string(jnode) + " at jnode " + std::to_string(global2local[uid]));
        }
    }

    std::vector<std::vector<int>> recv_halo( mpi_size );
    std::vector<std::vector<int>> recv_gidx( mpi_size );

    // send this data to the correct partition
    mpi::comm().allToAll( send_halo, recv_halo );
    mpi::comm().allToAll( send_gidx, recv_gidx );

    // adjust the values of the remote_halo field according to the results
    for ( idx_t p =0; p < mpi_size; ++p) {
      for ( idx_t i = 0; i < recv_halo[p].size(); ++i ) {
          //auto jnode = std::find(gidx.begin(), gidx.end(), recv_gidx[p][i]);
          idx_t found_idx = -1;
          gidx_t uid     = recv_gidx[p][i];
          Unique2Node::const_iterator found = global2local.find(uid);
          if (found != global2local.end()) {
              found_idx = found->second;
          }
          ATLAS_ASSERT(found_idx != -1,
              "global index not found: " + std::to_string(recv_gidx[p][i]));
          remote_halo (found_idx) = recv_halo[p][i];
      }
    }

}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
