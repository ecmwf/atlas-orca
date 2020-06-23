/*
 * (C) Copyright 2013 ECMWF.
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
#include "atlas/output/Gmsh.h"

#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test plugin" ) {
    // Because the library atlas_orca is not linked in, this test relies on dynamic
    // loading using environment variables:
    // ATLAS_PLUGINS_SEARCH_PATHS=<path-to-binary-dir>
    // ATLAS_PLUGINS=orca

    EXPECT( eckit::system::Library::exists( "atlas_orca" ) );
    EXPECT( bool( MeshGenerator( "orca" ) ) );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
