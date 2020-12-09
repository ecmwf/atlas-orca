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

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "atlas/util/Geometry.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

#include "atlas-orca/OrcaGrid.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

Grid cache_or_create( const std::string name ) {
    static std::map<std::string, Grid> grids;
    if ( grids.find( name ) == grids.end() ) {
        grids[name] = Grid( name );
    }
    return grids[name];
}

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

#if 1

size_t peakMemory() {
    return eckit::system::ResourceUsage().maxResidentSetSize();
}

CASE( "test generate orca mesh" ) {
    std::string gridname = "ORCA1_V";
    SECTION( gridname ) {
        auto grid = cache_or_create( gridname );
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

    SECTION( "auto_generate" ) { auto mesh = Mesh{cache_or_create( gridname )}; }
}

CASE( "test orca grid iterator" ) {
    struct Section {
        std::string gridname;
        size_t size;
    };

    std::vector<Section> sections{
        {"ORCA2_T", 27118},
        {"ORCA1_T", 105704},
        {"ORCA025_T", 1472282},
    };
    for ( auto& section : sections ) {
        std::string gridname = section.gridname;
        SECTION( gridname ) {
            OrcaGrid grid = cache_or_create( gridname );

            EXPECT_EQ( grid.size(), section.size );

            Log::info() << "grid.footprint() = " << eckit::Bytes( grid.footprint() ) << std::endl;

            idx_t n = 0;
            {
                auto trace = Trace( Here(), "iterating" );
                for ( auto& p : grid.lonlat() ) {
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
                EXPECT_EQ( mesh.nodes().size(), grid.size() );
            }

            // Now with extra virtual point at south pole
            EXPECT_EQ( ( MeshGenerator{"orca", util::Config( "include_pole", true )}.generate( grid ).nodes().size() ),
                       grid.size() + grid.nx() + grid->haloEast() + grid->haloWest() );
            EXPECT_EQ( ( MeshGenerator{"orca", util::Config( "force_include_south_pole", true )}
                             .generate( grid )
                             .nodes()
                             .size() ),
                       grid.size() + grid.nx() + grid->haloEast() + grid->haloWest() );
        }
    }
}
#endif


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


#if 1

CASE( "test orca mesh halo" ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",    //
        "ORCA1_T",    //
        "ORCA025_T",  //
    };
    auto& out = Log::debug();
    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            auto mesh = Mesh{cache_or_create( gridname )};
            REQUIRE( mesh.grid() );
            EXPECT( mesh.grid().name() == gridname );
            build_periodic_boundaries( mesh );
            auto remote_idx = array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() );
            auto ij         = array::make_view<idx_t, 2>( mesh.nodes().field( "ij" ) );
            idx_t count{0};

            functionspace::NodeColumns fs{mesh};
            Field field   = fs.createField<double>( option::name( "bla" ) );
            auto f        = array::make_view<double, 1>( field );
            OrcaGrid grid = mesh.grid();
            for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                if ( remote_idx( jnode ) < 0 ) {
                    auto p = PointIJ{ij( jnode, 0 ), ij( jnode, 1 )};
                    PointIJ master;
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
#endif


#if 1

CASE( "test periodicity " ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",     //
        "ORCA1_T",     //
        "eORCA1_T",    //
        "ORCA025_T",   //
        "eORCA025_T",  //
        // "ORCA12_T", //
        // "eORCA12_T",//
    };
    auto patch = std::map<std::string, bool>{
        {"ORCA2_T", false},     //
        {"ORCA1_T", true},      //
        {"eORCA1_T", true},     //
        {"ORCA025_T", false},   //
        {"eORCA025_T", false},  //
        {"ORCA12_T", false},    //
        {"eORCA12_T", false},   //
    };

    auto& out = Log::debug();
    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            OrcaGrid grid( cache_or_create( gridname ) );

            if ( true ) {
                idx_t jsubtract  = patch.at( gridname ) ? 1 : 2;
                double tolerance = 1.e-4;
                for ( idx_t j = grid.ny() - jsubtract; j < grid.ny() + grid->haloNorth(); ++j ) {
                    out << "----------------row-----------------" << std::endl;
                    bool pivot = false;
                    for ( idx_t i = -grid->haloWest(); i < grid.nx() + grid->haloEast(); ++i ) {
                        bool is_ghost = grid.ghost( i, j );
                        PointIJ folded;
                        grid->index2ij( grid->periodicIndex( i, j ), folded.i, folded.j );
                        bool is_pivot = folded.i == i && folded.j == j;
                        if ( not pivot ) {
                            if ( folded.i <= i && folded.i ) {
                                out << "===pivot===   ";
                                pivot = true;
                            }
                            if ( folded.i < i && folded.i ) {
                                out << std::endl;
                            }
                        }
                        out << PointIJ{i, j} << " --> " << folded << "     " << grid.lonlat( i, j ) << std::endl;
                        if ( is_pivot ) {
                            EXPECT( not is_ghost );
                        }
                        else if ( is_ghost ) {
                            EXPECT( not grid.ghost( folded.i, folded.j ) );
                        }
                        else {
                            EXPECT( grid.ghost( folded.i, folded.j ) );
                        }
                        EXPECT_APPROX_EQ( grid.lonlat( i, j ), grid.lonlat( folded.i, folded.j ), tolerance );
                    }
                }
            }


            Geometry geometry( 1. );
            auto print = [&]( idx_t j, idx_t i ) {
                Log::info() << grid->index( i, j ) + 1 << "\t" << std::string( grid.ghost( i, j ) ? "G" : " " ) << "\t"
                            << PointIJ{i, j} << " \t" << grid.lonlat( i, j ) << "\t"
                            << geometry.xyz( grid.lonlat( i, j ) ) << std::endl;
            };

            //            for( idx_t j=146; j<149; ++j ) {
            //              print(j,-1);
            //              print(j,0);
            //              print(j,1);
            //              print(j,178);
            //              print(j,179);
            //              print(j,180);
            //            }

            if ( false ) {
                auto ref_IJ = PointIJ{178, 147};
                Log::info() << "Looking for matches for " << ref_IJ << std::endl;
                auto ref   = grid.lonlat( ref_IJ.i, ref_IJ.j );
                bool match = false;
                for ( idx_t j = grid.ny() - 4; j < grid.ny() + grid->haloNorth(); ++j ) {
                    for ( idx_t i = -grid->haloWest(); i < grid.nx() + grid->haloEast(); ++i ) {
                        if ( grid.lonlat( i, j ) == ref || approx_eq( grid.lonlat( i, j ), ref ) ) {
                            Log::info() << "matches " << PointIJ{i, j} << std::endl;
                            match = true;
                            //Log::info() << "ix_pivot = " << (double(i+ref_IJ.i))/2. << std::endl;
                        }
                    }
                }
                EXPECT( match );
            }
        }
    }
}

#endif

/*
ORCA1_T
{290,-1} --> {289,0}     {73,50.0109}
[0] EXPECT_APPROX_EQ( grid.lonlat(i,j), grid.lonlat(folded.i,folded.j) ) FAILED @ test_orca.cc +440
[0]  --> lhs = {73.000000000000,50.010940551758}
[0]  --> rhs = {73.010848999023,50.010940551758}

{290,359} --> {289,0}     {73,50.0109}
[0] EXPECT_APPROX_EQ( grid.lonlat(i,j), grid.lonlat(folded.i,folded.j) ) FAILED @ test_orca.cc +440
[0]  --> lhs = {73.000000000000,50.010940551758}
[0]  --> rhs = {73.010848999023,50.010940551758}

*/

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
