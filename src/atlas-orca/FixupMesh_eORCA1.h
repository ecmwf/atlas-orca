
#pragma once


#include "atlas-orca/FixupMesh.h"

namespace atlas {
namespace meshgenerator {
namespace orca {

class Fixup_eORCA1 : public FixupMesh {
    void execute( MeshBuilder& ) const override;

protected:
    Fixup_eORCA1( const std::string name, const eckit::Parametrisation& config );

private:
    std::string gridname;
    bool include_south_pole_{false};
};

//----------------------------------------------------------------------------------------------------------------------

class Fixup_eORCA1_T : public Fixup_eORCA1 {
public:
    Fixup_eORCA1_T( const eckit::Parametrisation& config ) : Fixup_eORCA1( "eORCA1_T", config ) {}
};

class Fixup_eORCA1_F : public Fixup_eORCA1 {
public:
    Fixup_eORCA1_F( const eckit::Parametrisation& config ) : Fixup_eORCA1( "eORCA1_F", config ) {}
};

class Fixup_eORCA1_U : public Fixup_eORCA1 {
public:
    Fixup_eORCA1_U( const eckit::Parametrisation& config ) : Fixup_eORCA1( "eORCA1_U", config ) {}
};

class Fixup_eORCA1_V : public Fixup_eORCA1 {
public:
    Fixup_eORCA1_V( const eckit::Parametrisation& config ) : Fixup_eORCA1( "eORCA1_V", config ) {}
};

class Fixup_eORCA1_W : public Fixup_eORCA1 {
public:
    Fixup_eORCA1_W( const eckit::Parametrisation& config ) : Fixup_eORCA1( "eORCA1_W", config ) {}
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace orca
}  // namespace meshgenerator
}  // namespace atlas
