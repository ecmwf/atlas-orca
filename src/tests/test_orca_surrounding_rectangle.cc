/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "atlas-orca/meshgenerator/SurroundingRectangle.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/Point.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"

#include "tests/AtlasTestEnvironment.h"

using Grid = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

int wrap( idx_t value, idx_t lower, idx_t upper ) {
  // wrap around coordinate system when out of bounds
  const idx_t width = upper - lower;
  if (value < lower) {
    return wrap(value + width, lower, upper);
  }
  if (value > upper) {
    return wrap(value - width, lower, upper);
  }
  return value;
}

//-----------------------------------------------------------------------------

CASE("test surrounding rectangle ") {
  std::string gridname = "ORCA2_T";
  std::string distributionName = "checkerboard";

  auto rollup_plus = [](const double lon, const double lat) {
    return 1 + util::function::vortex_rollup(lon, lat, 0.0);
  };

  for (int halo = 0; halo < 3; ++halo) {
    SECTION(gridname + "_" + distributionName + "_halo" + std::to_string(halo)) {
      auto grid = OrcaGrid(gridname);
      auto partitioner_config = Config();
      partitioner_config.set("type", distributionName);
      auto partitioner = grid::Partitioner(partitioner_config);
      StructuredGrid::YSpace yspace{grid::LinearSpacing{
          {-80., 90.}, grid.ny(), true}};
      StructuredGrid::XSpace xspace{
          grid::LinearSpacing{{0., 360.}, grid.nx(), false}};
      StructuredGrid regular_grid{xspace, yspace};
      auto distribution = grid::Distribution(regular_grid, partitioner);

      orca::meshgenerator::SurroundingRectangle::Configuration cfg;
      cfg.mypart = mpi::rank();
      cfg.nparts = mpi::size();
      cfg.halosize = halo;
      cfg.nx_glb = grid.nx();
      cfg.ny_glb = grid.ny();
      std::cout << "[" << cfg.mypart << "] " << regular_grid.type() << std::endl;

      const idx_t cell_width = 1;
      if (regular_grid.ny() != grid.ny()) {
        std::cout << regular_grid.ny() << " != " << grid.ny() << std::endl;
      }
      EXPECT(regular_grid.ny() == grid.ny());
      for (idx_t ix = 0; ix < grid.ny(); ++ix) {
        if (regular_grid.nx(ix) != grid.nx()) {
          std::cout << regular_grid.nx(ix) << " != " << grid.nx() << std::endl;
        }
        EXPECT(regular_grid.nx(ix) == grid.nx());
      }

      orca::meshgenerator::SurroundingRectangle rectangle(distribution, cfg);

      auto regular_mesh = atlas::Mesh(regular_grid, partitioner);

      auto fview_lonlat =
          array::make_view<double, 2>(regular_mesh.nodes().lonlat());
      auto fview_glb_idx =
          array::make_view<gidx_t, 1>(regular_mesh.nodes().global_index());

      functionspace::NodeColumns regular_fs(regular_mesh);
      std::vector<int> indices;
      std::vector<bool> this_partition;
      for (uint64_t j = 0; j < rectangle.ny(); j++) {
        int iy_glb = rectangle.iy_min() + j;
        EXPECT(iy_glb < grid.ny() + cell_width + 2*halo);
        for (uint64_t i = 0; i < rectangle.nx(); i++) {
          int ix_glb = rectangle.ix_min() + i;
          EXPECT(ix_glb < grid.nx() + cell_width + 2*halo);
          auto ii = rectangle.index(i, j);
          indices.emplace_back(ii);
          atlas::orca::PointIJ ij = rectangle.global_periodic_ij(ix_glb, iy_glb);
          ATLAS_ASSERT(ij.i < cfg.nx_glb);
          ATLAS_ASSERT(ij.j < cfg.ny_glb);

          idx_t reg_grid_glb_idx  = regular_grid.index(ix_glb, iy_glb);
          idx_t orca_grid_glb_idx = grid.periodicIndex(ix_glb, iy_glb);
          idx_t reg_grid_remote_idx = 0;
          while(reg_grid_remote_idx < fview_glb_idx.size()) {
            if (reg_grid_glb_idx == fview_glb_idx(reg_grid_remote_idx))
              break;
            ++reg_grid_remote_idx;
          }

          if (rectangle.partition(ix_glb, iy_glb) == cfg.mypart) {
            this_partition.emplace_back(true);
          } else {
            this_partition.emplace_back(false);
          }
        }
      }

      int total_is_ghost =
          std::count(rectangle.is_ghost.begin(), rectangle.is_ghost.end(), true);

      EXPECT(indices.size() == rectangle.nx() * rectangle.ny());

      {
        // diagnostics
        auto total_on_partition =
            std::count(this_partition.begin(), this_partition.end(), true);
        auto not_on_partition =
            std::count(this_partition.begin(), this_partition.end(), false);

        //output::Gmsh gmsh(std::string("surroundingRect") +
        //                      std::to_string(cfg.nparts) + "_" + gridname + "_" +
        //                      distributionName + "_" + std::to_string(halo) +
        //                      ".msh",
        //                  Config("coordinates", "xy") | Config("info", true));
        //gmsh.write(regular_mesh);
      }

      std::cout << "(iy_min, iy_max) (" << rectangle.iy_min() << ", " << rectangle.iy_max() << ") "
                << "(ix_min, ix_max) (" << rectangle.ix_min() << ", " << rectangle.ix_max() << ") " << std::endl;
      if (cfg.nparts == 2) {
        if (cfg.mypart == 0) {
          EXPECT(rectangle.iy_min() == 0 - halo);
          EXPECT(rectangle.iy_max() == 147 + halo);
          EXPECT(rectangle.ix_min() == 0 - halo);
          EXPECT(rectangle.ix_max() == 90 + halo);
        }
        if (cfg.mypart == 1) {
          EXPECT(rectangle.iy_min() == 0 - halo);
          EXPECT(rectangle.iy_max() == 147 + halo);
          EXPECT(rectangle.ix_min() == 90 - halo);
          EXPECT(rectangle.ix_max() == 180 + halo);
        }
      }
      if (cfg.nparts == 4) {
        if (cfg.mypart == 0) {
          EXPECT(rectangle.iy_min() == 0 - halo);
          EXPECT(rectangle.iy_max() == 74 + halo);
          EXPECT(rectangle.ix_min() == 0 - halo);
          EXPECT(rectangle.ix_max() == 90 + halo);
        }
        if (cfg.mypart == 1) {
          EXPECT(rectangle.iy_min() == 0 - halo);
          EXPECT(rectangle.iy_max() == 74 + halo);
          EXPECT(rectangle.ix_min() == 89 - halo);
          EXPECT(rectangle.ix_max() == 180 + halo);
        }
        if (cfg.mypart == 2) {
          EXPECT(rectangle.iy_min() == 73 - halo);
          EXPECT(rectangle.iy_max() == 147 + halo);
          EXPECT(rectangle.ix_min() == 0 - halo);
          EXPECT(rectangle.ix_max() == 91 + halo);
        }
        if (cfg.mypart == 3) {
          EXPECT(rectangle.iy_min() == 73 - halo);
          EXPECT(rectangle.iy_max() == 147 + halo);
          EXPECT(rectangle.ix_min() == 90 - halo);
          EXPECT(rectangle.ix_max() == 180 + halo);
        }
      }
    }
  }
}

} // namespace test
} // namespace atlas

int main(int argc, char **argv) { return atlas::test::run(argc, argv); }
