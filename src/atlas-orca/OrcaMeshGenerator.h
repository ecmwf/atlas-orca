
#pragma once

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/util/Config.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}

namespace atlas {
class Mesh;
class OrcaGrid;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
}  // namespace grid
}  // namespace atlas
#endif

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class OrcaMeshGenerator : public MeshGenerator::Implementation {
public:
    OrcaMeshGenerator( const eckit::Parametrisation& = util::NoConfig() ) {}

    using MeshGenerator::Implementation::generate;

    void generate( const Grid&, const grid::Partitioner&, Mesh& ) const override;
    void generate( const Grid&, const grid::Distribution&, Mesh& ) const override;
    void generate( const Grid&, Mesh& ) const override;

    std::string type() const override { return "orca"; }
    static std::string static_type() { return "orca"; }

private:
    void hash( eckit::Hash& ) const override;

    util::Config options_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
