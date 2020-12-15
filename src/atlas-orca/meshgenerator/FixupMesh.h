
#pragma once

#include <memory>

#include "atlas/util/Config.h"

namespace atlas {
class Mesh;
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class MeshBuilder;

class FixupMesh {
public:
    static std::unique_ptr<FixupMesh> create( const eckit::Parametrisation& );
    void execute( Mesh& ) const;
    virtual void execute( MeshBuilder& ) const = 0;
    virtual ~FixupMesh() {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
