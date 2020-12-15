
#include "FixupMesh_ORCA.h"

#include "atlas/mesh.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

FixupMesh_ORCA::FixupMesh_ORCA( const eckit::Parametrisation& ) : FixupMesh_ORCA() {}

void FixupMesh_ORCA::execute( MeshBuilder& mesh_builder ) const {
    mesh_builder.mesh.metadata().set( "valid_elements", false );
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
