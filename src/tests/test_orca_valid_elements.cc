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
#include "atlas-orca/util/DetectInvalidElements.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;
using Quad3D = atlas::interpolation::element::Quad3D;


namespace atlas::test {

CASE( "test generate orca mesh" ) {
    std::vector<std::string> gridnames{
        "ORCA2_T", "ORCA2_F", "ORCA2_U", "ORCA2_V", "eORCA1_T", "eORCA1_F", "eORCA1_U", "eORCA1_V",
        "ORCA1_T", "ORCA1_F", "ORCA1_U", "ORCA1_V",
        "ORCA025_T", "ORCA025_F", "ORCA025_U", "ORCA025_V",
        "eORCA025_T", "eORCA025_F", "eORCA025_U", "eORCA025_V",
        "ORCA12_T", "ORCA12_F", "ORCA12_U", "ORCA12_V",
        "eORCA12_T", "eORCA12_F", "eORCA12_U", "eORCA12_V",
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
                                    << connectivity( e, 3 ) + 1 << " ] at lon/lats: " << pll[0] << " " << pll[1]
                                    << " " << pll[2] << " " << pll[3] << std::endl;
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

CASE( "test ORCA2 western Europe edge aspect" ) {
    const PointLonLat p_SW{ -1., 50. };
    const PointLonLat p_SE{ 1., 50. };
    const PointLonLat p_NE{ 1., 54. };
    const PointLonLat p_NW{ -1., 54. };

    orca::DetectInvalidElement::Statistics statistics;
    orca::DetectInvalidElement detect_orca2{ Config( "ORCA2", true ) };
    EXPECT( detect_orca2.invalid_element( p_SW, p_SE, p_NE, p_NW, statistics ) );
    EXPECT_EQUAL( statistics.western_europe_edge_aspect, 1 );

    orca::DetectInvalidElement detect_other_grid{ Config( "ORCA2", false ) };
    EXPECT( not detect_other_grid.invalid_element( p_SW, p_SE, p_NE, p_NW ) );

    orca::DetectInvalidElement detect_outside_region{ Config( "ORCA2", true ) };
    EXPECT( not detect_outside_region.invalid_element( PointLonLat{ 29., 50. }, PointLonLat{ 31., 50. },
                                                        PointLonLat{ 31., 54. }, PointLonLat{ 29., 54. } ) );
}

CASE( "test diagonal size is measured on the unit sphere" ) {
    orca::DetectInvalidElement detect{ util::NoConfig() };

    EXPECT( not detect.diagonal_too_large( PointLonLat{ 0., 89. }, PointLonLat{ 120., 89. },
                                           PointLonLat{ 120., 89.5 }, PointLonLat{ 0., 89.5 }, 5. ) );
    EXPECT( detect.diagonal_too_large( PointLonLat{ 0., 0. }, PointLonLat{ 10., 0. },
                                      PointLonLat{ 10., 10. }, PointLonLat{ 0., 10. }, 5. ) );
}

CASE( "test reported near-pole water elements which are self-intersecting from ORCA12_F" ) {
    struct Element {
        idx_t i;
        PointLonLat sw;
        PointLonLat se;
        PointLonLat ne;
        PointLonLat nw;
    };
    const std::array<Element, 8> elements{
        Element{ 1138, { 75.9516, 89.8581 }, { 71.9977, 89.8981 }, { 73.9155, 89.8981 }, { 70.0162, 89.8581 } },
        Element{ 1139, { 71.9977, 89.8981 }, { 75.4702, 89.9319 }, { 70.4494, 89.9319 }, { 73.9155, 89.8981 } },
        Element{ 1144, { -109.066, 89.9268 }, { -105.694, 89.893 }, { -107.525, 89.893 }, { -103.664, 89.9268 } },
        Element{ 1145, { -105.694, 89.893 }, { -109.753, 89.8533 }, { -104.02, 89.8533 }, { -107.525, 89.893 } },
        Element{ 3175, { -104.02, 89.8533 }, { -107.525, 89.893 }, { -105.694, 89.893 }, { -109.753, 89.8533 } },
        Element{ 3176, { -107.525, 89.893 }, { -103.664, 89.9268 }, { -109.066, 89.9268 }, { -105.694, 89.893 } },
        Element{ 3181, { 70.4494, 89.9319 }, { 73.9155, 89.8981 }, { 71.9977, 89.8981 }, { 75.4702, 89.9319 } },
        Element{ 3182, { 73.9155, 89.8981 }, { 70.0162, 89.8581 }, { 75.9516, 89.8581 }, { 71.9977, 89.8981 } },
    };

    orca::DetectInvalidElement detect{ util::NoConfig() };
    for ( const auto& element : elements ) {
        orca::DetectInvalidElement::Statistics statistics;
        EXPECT( not detect.invalid_quad_3d( element.sw, element.se, element.ne, element.nw ) );
        EXPECT( detect.invalid_element( element.sw, element.se, element.ne, element.nw, statistics ) );
        EXPECT_EQUAL( statistics.invalid_elements, 1 );
        EXPECT_EQUAL( statistics.self_intersecting, 1 );
        EXPECT( statistics.last_reason == orca::DetectInvalidElement::Reason::self_intersecting );
        EXPECT_EQUAL( statistics.diagonal_too_large, 0 );
        EXPECT_EQUAL( statistics.zero_diagonal, 0 );
        EXPECT_EQUAL( statistics.invalid_quads_3d, 0 );
    }
}


//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
