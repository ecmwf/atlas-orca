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

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/mesh/ElementType.h"

#include "atlas-orca/FixupMesh.h"
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

CASE( "test generate orca mesh" ) {
    using interpolation::element::Quad3D;

    bool modify          = true;
    std::string gridname = eckit::Resource<std::string>( "--grid", "ORCA2_T" );
    auto mesh            = Mesh{cache_or_create( gridname )};

    auto fixup_mesh =
        meshgenerator::FixupMesh::create( option::type( gridname ) | Config( "include_south_pole", true ) );
    if ( fixup_mesh ) {
        fixup_mesh->execute( mesh );
    }

    const auto& connectivity = mesh.cells().node_connectivity();
    auto lonlat              = array::make_view<double, 2>( mesh.nodes().lonlat() );
    const auto elem_glb_idx  = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
    auto flags               = array::make_view<int, 1>( mesh.cells().flags() );

    auto invalidate  = [&]( idx_t e ) { util::Topology::view( flags( e ) ).set( util::Topology::INVALID ); };
    auto invalidated = [&]( idx_t e ) { return util::Topology::view( flags( e ) ).check( util::Topology::INVALID ); };
    auto is_quad     = [&]( idx_t e ) { return mesh.cells().type_idx( e ) == 0; };

    Geometry geometry( 1. );
    bool has_invalid_quads = false;
    for ( idx_t e = 0; e < mesh.cells().size(); ++e ) {
        if ( not invalidated( e ) && is_quad( e ) ) {
            std::array<PointLonLat, 4> pll;
            std::array<PointXYZ, 4> pxyz;
            for ( idx_t n = 0; n < 4; ++n ) {
                pll[n]  = PointLonLat{lonlat( connectivity( e, n ), 0 ), lonlat( connectivity( e, n ), 1 )};
                pxyz[n] = geometry.xyz( pll[n] );
            }
            Quad3D quad{pxyz[0], pxyz[1], pxyz[2], pxyz[3]};
            if ( not quad.validate() ) {
                has_invalid_quads = true;
                //Log::info() << "Invalid quad [" << elem_glb_idx(e) << "] : [ " << pll[0] << ", " << pll[1] << ", " << pll[2] << ", " << pll[3] << " ]" << std::endl;
                Log::info() << "Invalid quad [" << elem_glb_idx( e ) << "] : [ " << connectivity( e, 0 ) + 1 << ", "
                            << connectivity( e, 1 ) + 1 << ", " << connectivity( e, 2 ) + 1 << ", "
                            << connectivity( e, 3 ) + 1 << " ]" << std::endl;
                //Log::info() << elem_glb_idx(e) << "," << std::endl;
                if ( modify ) {
                    invalidate( e );
                }
            }
        }
    }

    // Output mesh in different coordinates
    Config cfg;
    cfg.set( "info", true );
    cfg.set( "ghost", true );  //("water",true)("land",false);
    output::Gmsh( "ij.msh", Config( "coordinates", "ij" ) | cfg ).write( mesh );
    output::Gmsh( "lonlat.msh", Config( "coordinates", "lonlat" ) | cfg ).write( mesh );
    output::Gmsh( "xyz.msh", Config( "coordinates", "xyz" ) | cfg ).write( mesh );

    EXPECT( not has_invalid_quads );
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
