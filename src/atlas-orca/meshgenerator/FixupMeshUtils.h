#pragma once

#include <array>
#include <vector>

#include "atlas/array.h"
#include "atlas/mesh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

struct InvalidateElement {
    InvalidateElement( Mesh& mesh ) : flags( array::make_view<int, 1>( mesh.cells().flags() ) ) {}
    array::ArrayView<int, 1> flags;
    void operator()( gidx_t g ) { util::Topology::view( flags( g - 1 ) ).set( util::Topology::INVALID ); }
};

struct AddPoints {
    AddPoints( Mesh& m ) : mesh( m ), inode( mesh.nodes().size() ) {}
    idx_t operator()( const PointLonLat& p ) {
        points.push_back( p );
        return ++inode;
    }
    ~AddPoints() noexcept( false ) { add_points(); }
    void add_points() {
        if ( points.size() ) {
            ATLAS_DEBUG( "add_points" );
            idx_t n = mesh.nodes().size();
            mesh.nodes().resize( mesh.nodes().size() + points.size() );
            auto node_glb_idx = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );
            auto lonlat       = array::make_view<double, 2>( mesh.nodes().lonlat() );

            for ( const auto& p : points ) {
                lonlat( n, LON )  = p.lon();
                lonlat( n, LAT )  = p.lat();
                node_glb_idx( n ) = n + 1;
                ++n;
            }
            points.clear();
        }
    }
    Mesh mesh;
    idx_t inode;
    std::vector<PointLonLat> points;
};

struct AddElements {
    AddElements( Mesh& m ) : mesh( m ) {}
    struct Element {
        Element( gidx_t g1, gidx_t g2, gidx_t g3 ) {
            nodes[0] = g1 - 1;
            nodes[1] = g2 - 1;
            nodes[2] = g3 - 1;
        }
        Element( gidx_t g1, gidx_t g2, gidx_t g3, gidx_t g4 ) {
            nodes[0] = g1 - 1;
            nodes[1] = g2 - 1;
            nodes[2] = g3 - 1;
            nodes[3] = g4 - 1;
        }
        std::array<idx_t, 4> nodes;
        bool ghost = false;
    };
    Element& operator()( gidx_t g1, gidx_t g2, gidx_t g3 ) {
        triags.emplace_back( g1, g2, g3 );
        return triags.back();
    }
    Element& operator()( gidx_t g1, gidx_t g2, gidx_t g3, gidx_t g4 ) {
        quads.emplace_back( g1, g2, g3, g4 );
        return quads.back();
    }
    ~AddElements() noexcept( false ) { add_elements(); }

    void add_elements() {
        using util::Topology;
        if ( triags.empty() and quads.empty() ) {
            return;
        }
        idx_t e = mesh.cells().size();
        if ( triags.size() ) {
            Log::debug() << "Adding " << triags.size() << " triangles" << std::endl;
            mesh.cells().add( new mesh::temporary::Triangle(), triags.size() );
        }
        if ( quads.size() ) {
            mesh.cells().add( new mesh::temporary::Quadrilateral(), quads.size() );
            Log::debug() << "Adding " << quads.size() << " quads" << std::endl;
        }

        mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
        auto elem_halo                                        = array::make_view<int, 1>( mesh.cells().halo() );
        auto elem_glb_idx    = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
        auto elem_flags_view = array::make_view<int, 1>( mesh.cells().flags() );
        auto elem_flags      = [&]( idx_t i ) { return Topology::view( elem_flags_view( i ) ); };

        auto set_elem = [&]( idx_t e, const Element& element ) {
            node_connectivity.set( e, element.nodes.data() );
            elem_halo( e )    = 0;
            elem_glb_idx( e ) = e + 1;
            elem_flags( e ).reset();
            elem_flags( e ).set( Topology::LAND );
            if ( element.ghost ) {
                elem_flags( e ).set( Topology::GHOST );
            }
        };
        for ( const auto& element : triags ) {
            set_elem( e, element );
            ++e;
        }
        for ( const auto& element : quads ) {
            set_elem( e, element );
            ++e;
        }
        triags.clear();
        quads.clear();
    }
    Mesh mesh;
    std::vector<Element> triags;
    std::vector<Element> quads;
};

struct MeshBuilder {
    MeshBuilder( Mesh& m ) : mesh( m ), add_point{m}, add_element{m} {}
    Mesh mesh;
    AddPoints add_point;
    AddElements add_element;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
