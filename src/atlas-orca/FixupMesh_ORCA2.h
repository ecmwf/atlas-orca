
#pragma once


#include "atlas-orca/FixupMesh.h"

namespace atlas {
namespace meshgenerator {
namespace orca {

class Fixup_ORCA2 : public FixupMesh {
public:
    enum GridType { T,U,V,F };
    void execute( MeshBuilder& ) const override;

protected:
    Fixup_ORCA2( GridType grid_type, const eckit::Parametrisation& ) : type( grid_type ) {}

private:
    GridType type;
};

//----------------------------------------------------------------------------------------------------------------------

class Fixup_ORCA2_T : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_T( const eckit::Parametrisation& config ) : Fixup_ORCA2( T, config ) {}
};

class Fixup_ORCA2_F : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_F( const eckit::Parametrisation& config ) : Fixup_ORCA2( F, config ) {}
};

class Fixup_ORCA2_U : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_U( const eckit::Parametrisation& config ) : Fixup_ORCA2( U, config ) {}
};

class Fixup_ORCA2_V : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_V( const eckit::Parametrisation& config ) : Fixup_ORCA2( V, config ) {}
};

class Fixup_ORCA2_W : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_W( const eckit::Parametrisation& config ) : Fixup_ORCA2( T, config ) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace orca
}  // namespace meshgenerator
}  // namespace atlas
