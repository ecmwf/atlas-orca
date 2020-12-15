
#pragma once

#include <initializer_list>
#include <memory>

#include "atlas-orca/meshgenerator/FixupMesh.h"
#include "atlas-orca/meshgenerator/FixupMeshUtils.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class FixupMesh_ORCA : public FixupMesh {
public:
    enum GridType
    {
        T,
        F,
        U,
        V,
        W = T
    };

    FixupMesh_ORCA() = default;
    FixupMesh_ORCA( const eckit::Parametrisation& );
    void execute( MeshBuilder& ) const override;

protected:
    bool any_of( const GridType t, const std::initializer_list<GridType>& list ) const {
        for ( const auto& value : list ) {
            if ( t == value ) {
                return true;
            }
        }
        return false;
    }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
