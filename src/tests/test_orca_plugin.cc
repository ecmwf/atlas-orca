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
#include "atlas/meshgenerator.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas::test {

//-----------------------------------------------------------------------------

CASE( "test plugin" ) {
    // Because the library atlas_orca is not linked in, this test relies on dynamic
    // loading using environment variables:
    // ATLAS_PLUGINS_SEARCH_PATHS=<path-to-binary-dir>
    // ATLAS_PLUGINS=orca

    EXPECT( eckit::system::Library::exists( "atlas-orca" ) );
    EXPECT( bool( MeshGenerator( "orca" ) ) );
}

//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
