/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas::test {

namespace {

struct SplitCommunicator {
    SplitCommunicator() { mpi::comm( "world" ).split( mpi::comm( "world" ).rank() % 2, "split" ); }

    ~SplitCommunicator() {
        if ( eckit::mpi::hasComm( "split" ) ) {
            eckit::mpi::deleteComm( "split" );
        }
    }
};

}  // namespace

CASE( "ORCA mesh uses the partitioner communicator" ) {
    SplitCommunicator split_communicator;
    auto& split = mpi::comm( "split" );

    Grid grid{ mpi::comm( "world" ).rank() % 2 == 0 ? "ORCA2_T" : "eORCA1_T" };
    grid::Partitioner partitioner{
        "checkerboard", util::Config( "mpi_comm", "split" ) | util::Config( "partitions", split.size() ) };
    MeshGenerator meshgenerator{ "orca", util::Config( "halo", 1 ) };
    Mesh mesh = meshgenerator.generate( grid, partitioner );

    EXPECT_EQ( partitioner.mpi_comm(), "split" );
    EXPECT_EQ( mesh.mpi_comm(), "split" );
    EXPECT_EQ( mpi::comm().name(), "world" );
}

}  // namespace atlas::test

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}