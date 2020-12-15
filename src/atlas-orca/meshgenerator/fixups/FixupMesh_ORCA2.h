
#pragma once


#include "FixupMesh_ORCA.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

class FixupMesh_ORCA2 : public FixupMesh_ORCA {
public:
    void execute( MeshBuilder& ) const override;

protected:
    FixupMesh_ORCA2( GridType grid_type, const eckit::Parametrisation& ) : type( grid_type ) {}

private:
    GridType type;
};

//----------------------------------------------------------------------------------------------------------------------

class FixupMesh_ORCA2_T : public FixupMesh_ORCA2 {
public:
    FixupMesh_ORCA2_T( const eckit::Parametrisation& config ) : FixupMesh_ORCA2( T, config ) {}
};

class FixupMesh_ORCA2_F : public FixupMesh_ORCA2 {
public:
    FixupMesh_ORCA2_F( const eckit::Parametrisation& config ) : FixupMesh_ORCA2( F, config ) {}
};

class FixupMesh_ORCA2_U : public FixupMesh_ORCA2 {
public:
    FixupMesh_ORCA2_U( const eckit::Parametrisation& config ) : FixupMesh_ORCA2( U, config ) {}
};

class FixupMesh_ORCA2_V : public FixupMesh_ORCA2 {
public:
    FixupMesh_ORCA2_V( const eckit::Parametrisation& config ) : FixupMesh_ORCA2( V, config ) {}
};

class FixupMesh_ORCA2_W : public FixupMesh_ORCA2 {
public:
    FixupMesh_ORCA2_W( const eckit::Parametrisation& config ) : FixupMesh_ORCA2( W, config ) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
