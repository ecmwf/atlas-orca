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
#include "atlas/grid/Distribution.h"
#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/meshgenerator/SurroundingRectangle.h"


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}

namespace atlas {
class Mesh;
class OrcaGrid;
}  // namespace atlas


namespace atlas::grid {
class Distribution;
}  // namespace atlas::grid

#endif


namespace atlas::orca::meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class OrcaMeshGenerator : public MeshGenerator::Implementation {
public:
    explicit OrcaMeshGenerator( const eckit::Parametrisation& = util::NoConfig() );

    using MeshGenerator::Implementation::generate;

    void generate( const Grid&, const grid::Partitioner&, Mesh& ) const override;
    void generate( const Grid&, const grid::Distribution&, Mesh& ) const override;
    void generate( const Grid&, Mesh& ) const override;

    std::string type() const override { return "orca"; }
    static std::string static_type() { return "orca"; }

private:
    void hash( eckit::Hash& ) const override;
    void build_remote_index(Mesh& mesh) const;

    bool include_pole_{ false };
    bool fixup_{ true };
    int nparts_;
    int mypart_;
    int halosize_{0};
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas::orca::meshgenerator
