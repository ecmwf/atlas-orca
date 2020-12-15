
#pragma once


#include "FixupMesh_ORCA.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class FixupMesh_eORCA1 : public FixupMesh_ORCA {
public:
    void execute( MeshBuilder& ) const override;

protected:
    FixupMesh_eORCA1( GridType grid_type, const eckit::Parametrisation& config );

private:
    GridType type;
    bool include_south_pole_{false};
};

//----------------------------------------------------------------------------------------------------------------------

class FixupMesh_eORCA1_T : public FixupMesh_eORCA1 {
public:
    FixupMesh_eORCA1_T( const eckit::Parametrisation& config ) : FixupMesh_eORCA1( T, config ) {}
};

class FixupMesh_eORCA1_F : public FixupMesh_eORCA1 {
public:
    FixupMesh_eORCA1_F( const eckit::Parametrisation& config ) : FixupMesh_eORCA1( F, config ) {}
};

class FixupMesh_eORCA1_U : public FixupMesh_eORCA1 {
public:
    FixupMesh_eORCA1_U( const eckit::Parametrisation& config ) : FixupMesh_eORCA1( U, config ) {}
};

class FixupMesh_eORCA1_V : public FixupMesh_eORCA1 {
public:
    FixupMesh_eORCA1_V( const eckit::Parametrisation& config ) : FixupMesh_eORCA1( V, config ) {}
};

class FixupMesh_eORCA1_W : public FixupMesh_eORCA1 {
public:
    FixupMesh_eORCA1_W( const eckit::Parametrisation& config ) : FixupMesh_eORCA1( W, config ) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
