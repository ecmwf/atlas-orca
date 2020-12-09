
#pragma once

#include <memory>

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class MeshBuilder;

class FixupMesh {
public:
    static std::unique_ptr<FixupMesh> create( const util::Config& );
    void execute( Mesh& ) const;
    virtual void execute( MeshBuilder& ) const = 0;
    virtual ~FixupMesh() {}
};

//----------------------------------------------------------------------------------------------------------------------

namespace orca {

class FixupMesh_ORCA : public FixupMesh {
public:
    FixupMesh_ORCA( const eckit::Parametrisation& );
    void execute( MeshBuilder& ) const override;

private:
    bool include_south_pole_{false};
};

}  // namespace orca

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
