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

#include "eckit/exception/Exceptions.h" 
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Distribution.h"
#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}

namespace atlas {
class OrcaGrid;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
}  // namespace grid
}  // namespace atlas
#endif

namespace atlas::orca::meshgenerator {

//----------------------------------------------------------------------------------------------------------------------
//
class SurroundingRectangle {
 public:
    struct Configuration {
        int nparts;
        int mypart;
        int halosize;
        int nx_glb;
        int ny_glb;
        Configuration() :
            nparts(std::numeric_limits<int>::lowest())
            , mypart(std::numeric_limits<int>::lowest())
            , halosize(std::numeric_limits<int>::lowest())
            , nx_glb(std::numeric_limits<int>::lowest())
            , ny_glb(std::numeric_limits<int>::lowest()) {}
        void check_consistency() const {
            const auto check = [&](const int value) {
                if (value == std::numeric_limits<int>::lowest())
                  eckit::BadParameter("atlas-orca/meshgenerator/SurroundingRectangle: not all parameters set"); 
            };
            check(nparts);
            check(mypart);
            check(halosize);
            check(nx_glb);
            check(ny_glb);
        }
    };
    SurroundingRectangle(const grid::Distribution& distribution, const Configuration& cfg);
    PointIJ global_periodic_ij( idx_t ix, idx_t iy ) const;
    int index( int i, int j ) const;
    int partition( idx_t i, idx_t j ) const;
    int global_partition( idx_t ix_glb, idx_t iy_glb ) const;
    std::vector<int> parts;
    std::vector<int> halo;
    std::vector<int> is_ghost;
    uint64_t nx() const { return nx_; };
    uint64_t ny() const { return ny_; };
    int ix_min() const { return ix_min_; };
    int iy_min() const { return iy_min_; };
    int ix_max() const { return ix_max_; };
    int iy_max() const { return iy_max_; };
    int halosize() const { return halosize_; };
    uint64_t nb_real_nodes_owned_by_rectangle;

 private:
    const grid::Distribution distribution_;
    const OrcaGrid orca_;
    const Configuration cfg_;
    uint64_t nx_, ny_;
    int ix_min_, ix_max_;
    int iy_min_, iy_max_;
    int halosize_;
};
}  // namespace atlas::orca::meshgenerator 
