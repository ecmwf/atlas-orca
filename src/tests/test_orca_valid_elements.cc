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

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/mesh/ElementType.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;
using Quad3D = atlas::interpolation::element::Quad3D;


namespace atlas::test {

CASE( "test generate orca mesh" ) {
    std::vector<std::string> gridnames{
        "ORCA2_T", "ORCA2_F", "ORCA2_U", "ORCA2_V", "eORCA1_T", "eORCA1_F", "eORCA1_U", "eORCA1_V",
        "ORCA1_T", "ORCA1_F", "ORCA1_U", "ORCA1_V",
        "eORCA025_T", "eORCA025_F", "eORCA025_U", "eORCA025_V",
    };

    bool gmsh_output          = eckit::Resource<bool>( "--gmsh", false );
    std::string grid_resource = eckit::Resource<std::string>( "--grid", "" );
    if ( not grid_resource.empty() ) {
        gridnames = { grid_resource };
    }
    for ( const auto& gridname : gridnames ) {
        SECTION( gridname ) {
            auto mesh = Mesh{ gridname };

            const auto& connectivity = mesh.cells().node_connectivity();
            auto lonlat              = array::make_view<double, 2>( mesh.nodes().lonlat() );
            const auto elem_glb_idx  = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
            auto flags               = array::make_view<int, 1>( mesh.cells().flags() );

            auto invalidated = [&]( idx_t e ) {
                return util::Topology::view( flags( e ) ).check( util::Topology::INVALID );
            };

            geometry::Earth geometry;
            bool has_invalid_quads = false;
            for ( idx_t e = 0; e < mesh.cells().size(); ++e ) {
                if ( not invalidated( e ) ) {
                    std::array<PointLonLat, 4> pll;
                    std::array<PointXYZ, 4> pxyz;
                    for ( idx_t n = 0; n < 4; ++n ) {
                        pll[n]  = PointLonLat{ lonlat( connectivity( e, n ), 0 ), lonlat( connectivity( e, n ), 1 ) };
                        pxyz[n] = geometry.xyz( pll[n] );
                    }
                    Quad3D quad{ pxyz[0], pxyz[1], pxyz[2], pxyz[3] };
                    if ( not quad.validate() ) {
                        has_invalid_quads = true;
                        Log::info() << "Invalid quad [" << elem_glb_idx( e ) << "] : [ " << connectivity( e, 0 ) + 1
                                    << ", " << connectivity( e, 1 ) + 1 << ", " << connectivity( e, 2 ) + 1 << ", "
                                    << connectivity( e, 3 ) + 1 << " ]" << std::endl;
                    }
                }
            }

            if ( gmsh_output ) {
                // Output mesh in different coordinates
                Config cfg;
                cfg.set( "info", true );
                cfg.set( "ghost", true );  //("water",true)("land",false);
                output::Gmsh( gridname + "-ij.msh", Config( "coordinates", "ij" ) | cfg ).write( mesh );
                output::Gmsh( gridname + "-lonlat.msh", Config( "coordinates", "lonlat" ) | cfg ).write( mesh );
                output::Gmsh( gridname + "-xyz.msh", Config( "coordinates", "xyz" ) | cfg ).write( mesh );
            }
            EXPECT( not has_invalid_quads );
        }
    }
}


//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
