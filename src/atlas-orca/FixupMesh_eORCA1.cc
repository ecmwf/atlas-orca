
#include "FixupMesh_eORCA1.h"

#include "FixupMeshUtils.h"

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------


namespace {
using GridType = orca::Fixup_eORCA1::GridType;

bool any_of( const GridType t, const std::initializer_list<GridType>& list ) {
    for ( const auto& value : list ) {
        if ( t == value ) {
            return true;
        }
    }
    return false;
}

}  // namespace

void orca::Fixup_eORCA1::execute( MeshBuilder& mesh_builder ) const {
    ATLAS_TRACE( "FixupMesh " + type );

    auto mesh = mesh_builder.mesh;

    bool valid_elements{false};
    mesh.metadata().get( "valid_elements", valid_elements );
    if ( valid_elements ) {
        // Nothing to be done here
        return;
    }
    mesh.metadata().set( "valid_elements", true );

    bool includes_south_pole{false};
    mesh.metadata().get( "includes_south_pole", includes_south_pole );

    InvalidateElement invalidate( mesh );
    auto& add_element = mesh_builder.add_element;

    for ( gidx_t g = 0; g < 48; ++g ) {
        invalidate( 200 + g * 361 );
    }

    if ( any_of( type, {T} ) ) {
        for ( gidx_t g = 0; g < 41; ++g ) {
            invalidate( 40 + g * 361 );
        }
        invalidate( 14841 );
        invalidate( 15202 );
        invalidate( 15563 );
        invalidate( 15924 );
        invalidate( 16285 );
        invalidate( 16646 );
        invalidate( 17528 );
    }
    if ( any_of( type, {U, F} ) ) {
        for ( gidx_t g = 0; g < 45; ++g ) {
            invalidate( 39 + g * 361 );
        }
    }
    if ( any_of( type, {V} ) ) {
        for ( gidx_t g = 0; g < 45; ++g ) {
            invalidate( 40 + g * 361 );
        }
    }
    if ( any_of( type, {T} ) ) {
        add_element( 17577, 17939, 17938 );
        for ( idx_t i = 0; i > ( -20 ); --i ) {
            add_element( 17055 + 362 * ( i - 1 ), 17055 + 362 * i, 17054 );
        }
        add_element( 17054, 16692, 9815 );
        for ( idx_t i = 0; i > ( -6 ); --i ) {
            add_element( 9815 + 362 * ( i - 1 ), 9815 + 362 * i, 16692 );
        }
        add_element( 16692, 16330, 7643 );

        if ( includes_south_pole ) {
            // left fold
            invalidate( 119531 );
            invalidate( 119532 );
            invalidate( 119533 );

            // right fold
            for ( gidx_t g = 119658; g <= 119691; ++g ) {
                invalidate( g );
            }
        }
    }
    if ( any_of( type, {U} ) ) {
        for ( int i = 0; i < 26; ++i ) {
            add_element( 16329, 6556 + i * 362, 6918 + ( i + 1 ) * 362 );
        }
        add_element( 17214, 17577, 17576 );

        if ( includes_south_pole ) {
            Log::warning() << "Todo: add elements " << Here() << std::endl;
            mesh.metadata().set( "valid_elements", false );

        }
    }
    if ( any_of( type, {F} ) ) {
        for ( int i = 0; i < 26; ++i ) {
            add_element( 16329, 6556 + i * 362, 6918 + ( i + 1 ) * 362 );
        }
        add_element( 17214, 17577, 17576 );

        if ( includes_south_pole ) {
            Log::warning() << "Todo: add elements " << Here() << std::endl;
            mesh.metadata().set( "valid_elements", false );


        }
    }
    if ( any_of( type, {V} ) ) {
        for ( int i = 0; i < 26; ++i ) {
            add_element( 16330, 6919 + i * 362, 7281 + i * 362 );
        }
        add_element( 17215, 17577, 17576 );
        if ( includes_south_pole ) {
            Log::warning() << "Todo: add elements " << Here() << std::endl;
            mesh.metadata().set( "valid_elements", false );


        }
    }
}

orca::Fixup_eORCA1::Fixup_eORCA1( GridType grid_type, const eckit::Parametrisation& config ) : type(grid_type) {
    config.get( "include_south_pole", include_south_pole_ );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
