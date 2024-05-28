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

#include <limits>
#include <memory>

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Distribution.h"
#include "atlas/util/Point.h"
#include "atlas/util/Bitflags.h"
#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/meshgenerator/SurroundingRectangle.h"
#include "atlas-orca/util/PointIJ.h"


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}
#endif

namespace atlas::orca::meshgenerator {

//----------------------------------------------------------------------------------------------------------------------
class LocalOrcaGrid {
 public:
    std::vector<int> parts;
    std::vector<int> halo;
    std::vector<int> is_ghost_including_orca_halo;
    std::vector<int> is_ghost;
    std::vector<int> is_node;
    std::vector<int> is_cell;
    uint64_t size() const {return size_;}
    int ix_min() const {return ix_orca_min_;}
    int ix_max() const {return ix_orca_max_;}
    int iy_min() const {return iy_orca_min_;}
    int iy_max() const {return iy_orca_max_;}
    uint64_t nx() const {return nx_orca_;}
    uint64_t ny() const {return ny_orca_;}
    // number of real nodes on this partition
    uint64_t nb_used_real_nodes() const {return nb_used_real_nodes_;}
    // number of ghost nodes on this partition
    uint64_t nb_used_ghost_nodes() const {return nb_used_ghost_nodes_;}
    // number of nodes used by cells on this partition
    uint64_t nb_used_nodes() const {return nb_used_nodes_;}
    // number of cells on this partition
    uint64_t nb_cells() const {return nb_cells_;}
    // number of ghost cells on this partition
    uint64_t nb_used_ghost_cells() const {return nb_used_ghost_cells_;}
    // number of real cells on this paritition
    uint64_t nb_used_real_cells() const {return nb_used_real_cells_;}

    int index( idx_t ix, idx_t iy ) const;
    LocalOrcaGrid( const OrcaGrid& grid, const SurroundingRectangle& rectangle );
    PointIJ global_ij( idx_t ix, idx_t iy ) const;
    const PointXY grid_xy( idx_t ix, idx_t iy ) const;
    PointXY normalised_grid_xy( idx_t ix, idx_t iy ) const;
    PointIJ master_global_ij( idx_t ix, idx_t iy ) const;
    gidx_t master_global_index( idx_t ix, idx_t iy ) const;
    PointLonLat normalised_grid_master_lonlat( idx_t ix, idx_t iy ) const;
    PointIJ orca_haloed_global_grid_ij( idx_t ix, idx_t iy ) const;
    idx_t orca_haloed_global_grid_index( idx_t ix, idx_t iy ) const;
    void flags( idx_t ix, idx_t iy, util::detail::BitflagsView<int>& flag_view ) const;
    bool water( idx_t ix, idx_t iy ) const;
    bool orca_halo( idx_t ix, idx_t iy ) const;

 private:
    const OrcaGrid orca_;
    uint64_t size_;
    int ix_orca_min_;
    int ix_orca_max_;
    int iy_orca_min_;
    int iy_orca_max_;
    uint64_t nx_orca_;
    uint64_t ny_orca_;
    uint64_t nb_used_nodes_;
    uint64_t nb_used_real_nodes_;
    uint64_t nb_used_ghost_nodes_;
    uint64_t nb_cells_;
    uint64_t nb_used_ghost_cells_;
    uint64_t nb_used_real_cells_;
    double lon00_;
    util::NormaliseLongitude lon00_normaliser_;
};
}  // namespace atlas::orca::meshgenerator
