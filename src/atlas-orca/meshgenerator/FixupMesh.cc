
#include "FixupMesh.h"
#include "FixupMeshUtils.h"

#include "fixups/FixupMesh_ORCA2.h"
#include "fixups/FixupMesh_eORCA1.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas/mesh.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

namespace {
template <class T, class... Args>
std::unique_ptr<T> make_unique( Args&&... args ) {
    return std::unique_ptr<T>( new T( std::forward<Args>( args )... ) );
}
}  // namespace

std::unique_ptr<FixupMesh> FixupMesh::create( const eckit::Parametrisation& config ) {
    std::string type{};
    config.get( "type", type );
    if ( type == "ORCA2_T" ) {
        return make_unique<FixupMesh_ORCA2_T>( config );
    }
    if ( type == "ORCA2_F" ) {
        return make_unique<FixupMesh_ORCA2_F>( config );
    }
    if ( type == "ORCA2_U" ) {
        return make_unique<FixupMesh_ORCA2_U>( config );
    }
    if ( type == "ORCA2_V" ) {
        return make_unique<FixupMesh_ORCA2_V>( config );
    }
    if ( type == "ORCA2_W" ) {
        return make_unique<FixupMesh_ORCA2_W>( config );
    }
    if ( type == "eORCA1_T" ) {
        return make_unique<FixupMesh_eORCA1_T>( config );
    }
    if ( type == "eORCA1_F" ) {
        return make_unique<FixupMesh_eORCA1_F>( config );
    }
    if ( type == "eORCA1_U" ) {
        return make_unique<FixupMesh_eORCA1_U>( config );
    }
    if ( type == "eORCA1_V" ) {
        return make_unique<FixupMesh_eORCA1_V>( config );
    }
    if ( type == "eORCA1_W" ) {
        return make_unique<FixupMesh_eORCA1_W>( config );
    }
    return make_unique<FixupMesh_ORCA>( config );
}

void FixupMesh::execute( Mesh& mesh ) const {
    MeshBuilder mesh_builder( mesh );
    execute( mesh_builder );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
