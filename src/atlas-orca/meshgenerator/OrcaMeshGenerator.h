/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

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
namespace orca {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class OrcaMeshGenerator : public MeshGenerator::Implementation {
public:
    OrcaMeshGenerator( const eckit::Parametrisation& = util::NoConfig() );

    using MeshGenerator::Implementation::generate;

    void generate( const Grid&, const grid::Partitioner&, Mesh& ) const override;
    void generate( const Grid&, const grid::Distribution&, Mesh& ) const override;
    void generate( const Grid&, Mesh& ) const override;

    std::string type() const override { return "orca"; }
    static std::string static_type() { return "orca"; }

private:
    void hash( eckit::Hash& ) const override;

    bool include_pole_{false};
    bool fixup_{true};
    int halosize_{0};
    int nparts_;
    int mypart_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
