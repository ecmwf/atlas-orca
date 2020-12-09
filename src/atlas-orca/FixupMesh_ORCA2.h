
#pragma once


#include "atlas-orca/FixupMesh.h"

namespace atlas {
namespace meshgenerator {
namespace orca {

class Fixup_ORCA2 : public FixupMesh {
    void execute( MeshBuilder& ) const override;

protected:
    Fixup_ORCA2( const std::string name, const eckit::Parametrisation& ) : gridname( name ) {}

private:
    std::string gridname;
};

//----------------------------------------------------------------------------------------------------------------------

class Fixup_ORCA2_T : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_T( const eckit::Parametrisation& config ) : Fixup_ORCA2( "ORCA2_T", config ) {}
};

class Fixup_ORCA2_F : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_F( const eckit::Parametrisation& config ) : Fixup_ORCA2( "ORCA2_F", config ) {}
};

class Fixup_ORCA2_U : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_U( const eckit::Parametrisation& config ) : Fixup_ORCA2( "ORCA2_U", config ) {}
};

class Fixup_ORCA2_V : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_V( const eckit::Parametrisation& config ) : Fixup_ORCA2( "ORCA2_V", config ) {}
};

class Fixup_ORCA2_W : public Fixup_ORCA2 {
public:
    Fixup_ORCA2_W( const eckit::Parametrisation& config ) : Fixup_ORCA2( "ORCA2_T", config ) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace orca
}  // namespace meshgenerator
}  // namespace atlas
