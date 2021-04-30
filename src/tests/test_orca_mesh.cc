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

void build_periodic_boundaries( Mesh& mesh );  // definition below

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
            test::build_periodic_boundaries( mesh );
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

// WORK IN PROGRESS: The routine in Atlas should be incorporating the following when
// use in parallel is required.

void build_periodic_boundaries( Mesh& mesh ) {
    using util::LonLatMicroDeg;
    using util::PeriodicTransform;
    using util::Topology;

    ATLAS_TRACE();
    bool periodic = false;
    mesh.metadata().get( "periodic", periodic );

    auto mpi_size = mpi::size();
    auto mypart   = mpi::rank();

    if ( !periodic ) {
        mesh::Nodes& nodes = mesh.nodes();

        auto flags_v = array::make_view<int, 1>( nodes.flags() );
        auto ridx    = array::make_indexview<idx_t, 1>( nodes.remote_index() );
        auto part    = array::make_view<int, 1>( nodes.partition() );
        auto ghost   = array::make_view<int, 1>( nodes.ghost() );

        int nb_nodes = nodes.size();

        auto xy = array::make_view<double, 2>( nodes.xy() );
        auto ij = array::make_view<idx_t, 2>( nodes.field( "ij" ) );

        // Identify my master and slave nodes on own partition
        // master nodes are at x=0,  slave nodes are at x=2pi
        std::map<uid_t, int> master_lookup;
        std::map<uid_t, int> slave_lookup;
        std::vector<int> master_nodes;
        master_nodes.reserve( 3 * nb_nodes );
        std::vector<int> slave_nodes;
        slave_nodes.reserve( 3 * nb_nodes );

        auto collect_slave_and_master_nodes = [&]() {
            for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
                auto flags = Topology::view( flags_v( jnode ) );

                if ( flags.check_all( Topology::BC | Topology::WEST ) ) {
                    flags.set( Topology::PERIODIC );
                    if ( part( jnode ) == mypart ) {
                        LonLatMicroDeg ll( xy( jnode, XX ), xy( jnode, YY ) );
                        master_lookup[ll.unique()] = jnode;
                        master_nodes.push_back( ll.lon() );
                        master_nodes.push_back( ll.lat() );
                        master_nodes.push_back( jnode );
                    }
                    // Log::info() << "master " << jnode << "  " << PointXY{ij( jnode, 0 ), ij( jnode, 1 )} << std::endl;
                }
                else if ( flags.check( Topology::BC | Topology::EAST ) ) {
                    flags.set( Topology::PERIODIC | Topology::GHOST );
                    ghost( jnode ) = 1;
                    LonLatMicroDeg ll( xy( jnode, XX ), xy( jnode, YY ) );
                    slave_lookup[ll.unique()] = jnode;
                    slave_nodes.push_back( ll.lon() );
                    slave_nodes.push_back( ll.lat() );
                    slave_nodes.push_back( jnode );
                    ridx( jnode ) = -1;
                    // Log::info() << "slave " << jnode << "  " << PointXY{ij( jnode, 0 ), ij( jnode, 1 )} << std::endl;
                }
            }
        };

        auto collect_slave_and_master_nodes_halo = [&]() {
            auto master_glb_idx = array::make_view<gidx_t, 1>( mesh.nodes().field( "master_global_index" ) );
            auto glb_idx        = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );
            for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
                auto flags = Topology::view( flags_v( jnode ) );
                if ( flags.check( Topology::PERIODIC ) ) {
                    if ( master_glb_idx( jnode ) == glb_idx( jnode ) ) {
                        if ( part( jnode ) == mypart ) {
                            master_lookup[master_glb_idx( jnode )] = jnode;
                            master_nodes.push_back( master_glb_idx( jnode ) );
                            master_nodes.push_back( 0 );
                            master_nodes.push_back( jnode );
                            //Log::info() << "master " << jnode << "  " << PointXY{ ij(jnode,0), ij(jnode,1)} << std::endl;
                        }
                    }
                    else {
                        slave_lookup[master_glb_idx( jnode )] = jnode;
                        slave_nodes.push_back( master_glb_idx( jnode ) );
                        slave_nodes.push_back( 0 );
                        slave_nodes.push_back( jnode );
                        ridx( jnode ) = -1;
                        //Log::info() << "slave " << jnode << "  " << PointXY{ ij(jnode,0), ij(jnode,1)} << std::endl;
                    }
                }
            }
        };

        bool require_PeriodicTransform;
        if ( mesh.nodes().has_field( "master_global_index" ) ) {
            collect_slave_and_master_nodes_halo();
            require_PeriodicTransform = false;
        }
        else {
            collect_slave_and_master_nodes();
            require_PeriodicTransform = true;
        }
        std::vector<std::vector<int>> found_master( mpi_size );
        std::vector<std::vector<int>> send_slave_idx( mpi_size );

        // Find masters on other tasks to send to me
        {
            int sendcnt = slave_nodes.size();
            std::vector<int> recvcounts( mpi_size );

            ATLAS_TRACE_MPI( ALLGATHER ) { mpi::comm().allGather( sendcnt, recvcounts.begin(), recvcounts.end() ); }

            std::vector<int> recvdispls( mpi_size );
            recvdispls[0] = 0;
            int recvcnt   = recvcounts[0];
            for ( idx_t jproc = 1; jproc < mpi_size; ++jproc ) {
                recvdispls[jproc] = recvdispls[jproc - 1] + recvcounts[jproc - 1];
                recvcnt += recvcounts[jproc];
            }
            std::vector<int> recvbuf( recvcnt );

            ATLAS_TRACE_MPI( ALLGATHER ) {
                mpi::comm().allGatherv( slave_nodes.begin(), slave_nodes.end(), recvbuf.begin(), recvcounts.data(),
                                        recvdispls.data() );
            }

            PeriodicTransform transform;
            for ( idx_t jproc = 0; jproc < mpi_size; ++jproc ) {
                found_master.reserve( master_nodes.size() );
                send_slave_idx.reserve( master_nodes.size() );
                array::LocalView<int, 2> recv_slave( recvbuf.data() + recvdispls[jproc],
                                                     array::make_shape( recvcounts[jproc] / 3, 3 ) );
                for ( idx_t jnode = 0; jnode < recv_slave.shape( 0 ); ++jnode ) {
                    uid_t slave_uid;
                    if ( require_PeriodicTransform ) {
                        LonLatMicroDeg slave( recv_slave( jnode, LON ), recv_slave( jnode, LAT ) );
                        transform( slave, -1 );
                        uid_t slave_uid = slave.unique();
                    }
                    else {
                        slave_uid = recv_slave( jnode, 0 );
                    }
                    //int debug_slave_idx  = recv_slave( jnode, 2 );
                    //ATLAS_DEBUG("slave " <<  (PointXY{ ij(debug_slave_idx,0), ij(debug_slave_idx,1)}) << "  uid = " << slave_uid);
                    if ( master_lookup.count( slave_uid ) ) {
                        int master_idx = master_lookup[slave_uid];
                        int slave_idx  = recv_slave( jnode, 2 );
                        //ATLAS_DEBUG( "  found in master  " << (PointXY{ij(master_idx,0),ij(master_idx,1)}) );
                        found_master[jproc].push_back( master_idx );
                        send_slave_idx[jproc].push_back( slave_idx );
                    }
                }
            }
        }

        // Fill in data to communicate
        std::vector<std::vector<int>> recv_slave_idx( mpi_size );
        std::vector<std::vector<int>> send_master_part( mpi_size );
        std::vector<std::vector<int>> recv_master_part( mpi_size );
        std::vector<std::vector<int>> send_master_ridx( mpi_size );
        std::vector<std::vector<int>> recv_master_ridx( mpi_size );

        //  std::vector< std::vector<int> > send_slave_part( mpi_size );
        //  std::vector< std::vector<int> > recv_slave_part( mpi_size );
        //  std::vector< std::vector<int> > send_slave_ridx( mpi_size );
        //  std::vector< std::vector<int> > recv_slave_ridx( mpi_size );

        {
            for ( idx_t jproc = 0; jproc < mpi_size; ++jproc ) {
                idx_t nb_found_master = static_cast<idx_t>( found_master[jproc].size() );
                send_master_part[jproc].resize( nb_found_master );
                send_master_ridx[jproc].resize( nb_found_master );
                for ( idx_t jnode = 0; jnode < nb_found_master; ++jnode ) {
                    int loc_idx                    = found_master[jproc][jnode];
                    send_master_part[jproc][jnode] = part( loc_idx );
                    send_master_ridx[jproc][jnode] = loc_idx;
                }

                //      int nb_found_slaves = found_slave[jproc].size();
                //      send_slave_glb_idx[jproc].resize(nb_found_slaves);
                //      send_slave_part   [jproc].resize(nb_found_slaves);
                //      send_slave_ridx   [jproc].resize(nb_found_slaves);
                //      for( int jnode=0; jnode<nb_found_slaves; ++jnode )
                //      {
                //        int loc_idx = found_slave[jproc][jnode];
                //        send_slave_glb_idx[jproc][jnode] = glb_idx( loc_idx );
                //        send_slave_part   [jproc][jnode] = part   ( loc_idx );
                //        send_slave_ridx   [jproc][jnode] = loc_idx;
                //      }
            }
        }

        // Communicate
        ATLAS_TRACE_MPI( ALLTOALL ) {
            mpi::comm().allToAll( send_slave_idx, recv_slave_idx );
            mpi::comm().allToAll( send_master_part, recv_master_part );
            mpi::comm().allToAll( send_master_ridx, recv_master_ridx );
            //  mpi::comm().allToAll( send_slave_part,     recv_slave_part    );
            //  mpi::comm().allToAll( send_slave_loc,      recv_slave_ridx    );
        }

        // Fill in periodic
        // unused // int nb_recv_master = 0;
        for ( idx_t jproc = 0; jproc < mpi_size; ++jproc ) {
            idx_t nb_recv = static_cast<idx_t>( recv_slave_idx[jproc].size() );
            for ( idx_t jnode = 0; jnode < nb_recv; ++jnode ) {
                idx_t slave_idx   = recv_slave_idx[jproc][jnode];
                part( slave_idx ) = recv_master_part[jproc][jnode];
                ridx( slave_idx ) = recv_master_ridx[jproc][jnode];
            }
        }
    }
    mesh.metadata().set( "periodic", true );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
