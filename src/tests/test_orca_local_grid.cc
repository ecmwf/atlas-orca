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
#include "atlas-orca/meshgenerator/LocalOrcaGrid.h"
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
#include "atlas/util/Bitflags.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"

#include "tests/AtlasTestEnvironment.h"

using Grid = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test surrounding local_orca ") {
  std::string gridname = "ORCA2_T";
  std::string distributionName = "checkerboard";

  auto rollup_plus = [](const double lon, const double lat) {
    return 1 + util::function::vortex_rollup(lon, lat, 0.0);
  };

  for (int halo = 0; halo <= 2; ++halo) {
    if ( ( (mpi::size() == 1) && (halo > 0) ) ||
         ( (mpi::size() == 1) && (distributionName != "serial") ) )
      continue;
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
      orca::meshgenerator::SurroundingRectangle rectangle(distribution, cfg);
      std::cout << "[" << cfg.mypart << "] rectangle.ix_min " <<  rectangle.ix_min()
                << " rectangle.ix_max " <<  rectangle.ix_max()
                << " rectangle.iy_min " <<  rectangle.iy_min()
                << " rectangle.iy_max " <<  rectangle.iy_max() << std::endl;
      orca::meshgenerator::LocalOrcaGrid local_orca(grid, rectangle);

      std::vector<idx_t> indices;
      std::vector<bool> this_partition;
      int inode_ghost = local_orca.nb_used_real_nodes();  // orca ghost nodes start counting after nonghost nodes
      int inode_nonghost = 0;
      for (uint64_t j = 0; j < local_orca.ny(); j++) {
        int iy_glb = local_orca.iy_min() + j;
        ASSERT_MSG(iy_glb < (grid.ny() + std::max(grid.haloNorth(), halo + 1)), "iy_glb in range " +
                   std::to_string(j) + " " + std::to_string(local_orca.iy_min()) + " " +
                   std::to_string(iy_glb) + " " + std::to_string(grid.ny()) + " " +
                   std::to_string(grid.haloNorth()) + " " + std::to_string(halo));
        for (uint64_t i = 0; i < local_orca.nx(); i++) {
          int ix_glb = local_orca.ix_min() + i;
          ASSERT_MSG(ix_glb < (grid.nx() + std::max(grid.haloEast(), halo + 1)), "ix_glb in range " +
                     std::to_string(i) + " " + std::to_string(local_orca.ix_min()) + " " +
                     std::to_string(ix_glb) + " " + std::to_string(grid.nx()) + " " +
                     std::to_string(grid.haloEast()) + " " + std::to_string(halo));
          auto ii = local_orca.index(i, j);
          indices.emplace_back(ii);

          gidx_t ohgg_idx = local_orca.orca_haloed_global_grid_index( i, j );
          const auto ohgg_ij = local_orca.orca_haloed_global_grid_ij( i, j );
          ASSERT_MSG(ohgg_ij.i < grid.nx() + grid.haloEast(),
             std::string("ohgg_ij.i ") + std::to_string(ohgg_ij.i)
             + std::string(" grid.nx()) ") + std::to_string(grid.nx())
             + std::string(" grid.haloEast() ") + std::to_string(grid.haloEast()) );
          ASSERT_MSG(ohgg_ij.j < grid.ny() + grid.haloNorth(),
             std::string("ohgg_ij.j ") + std::to_string(ohgg_ij.j)
             + std::string(" grid.ny() ") + std::to_string(grid.ny())
             + std::string(" grid.haloNorth() ") + std::to_string(grid.haloNorth()) );
          ASSERT_MSG(ohgg_ij.i >= -grid.haloWest(),
             std::string("ohgg_ij.i ") + std::to_string(ohgg_ij.i)
             + std::string(" grid.haloWest() ") + std::to_string(grid.haloWest()) );
          ASSERT_MSG(ohgg_ij.j >= -grid.haloSouth(),
             std::string("ohgg_ij.j ") + std::to_string(ohgg_ij.j)
             + std::string(" grid.haloSouth() ") + std::to_string(grid.haloSouth()) );

          const auto grid_xy        = local_orca.grid_xy( i, j );
          const auto normed_grid_xy = local_orca.normalised_grid_xy( i, j );

          if (halo == 0) {
            const auto ij_glb = local_orca.global_ij( i, j );
            const auto ij_glb_haloed = local_orca.orca_haloed_global_grid_ij( i, j );
            ASSERT_MSG(ij_glb.i == ij_glb_haloed.i,
               std::string("ij_glb.i != ij_glb_haloed.i ") + std::to_string(ij_glb.i)
               + std::string(" != ") + std::to_string(ij_glb_haloed.i));
            ASSERT_MSG(ij_glb.j == ij_glb_haloed.j,
               std::string("ij_glb.j != ij_glb_haloed.j ") + std::to_string(ij_glb.j)
               + std::string(" != ") + std::to_string(ij_glb_haloed.j));
          }

          if (local_orca.parts.at(ii) == cfg.mypart) {
            this_partition.emplace_back(true);
          } else {
            this_partition.emplace_back(false);
          }
          const auto water = local_orca.water(i, j);
          const auto halo = local_orca.halo.at(ii);
          int flags = 0;
          util::detail::BitflagsView<int> flags_view(flags);
          local_orca.flags( i, j, flags_view );

          // check points in the orca halo behave as expected.
          if (local_orca.orca_edge( i, j ) &&
              ((ix_glb >= grid.nx()) || (ix_glb > grid.nx()/2 && iy_glb > grid.ny()))) {
            // this grid point should be a ghost point.
            if ( ohgg_ij.j >= 0 or ohgg_ij.i < 0 or ohgg_ij.i >= grid.nx() ) {
              EXPECT(local_orca.is_ghost_including_orca_halo.at(ii)
                     >= grid.ghost( ohgg_ij.i, ohgg_ij.j ));
            }
            if ( (ix_glb >= -grid.haloWest()) && (ix_glb < grid.nx() + grid.haloEast()) &&
                 (iy_glb >= -grid.haloSouth()) && (iy_glb < grid.ny() + grid.haloNorth())) {
              // this grid point should not be a master grid point.
              EXPECT(local_orca.master_global_index(i, j)
                     != local_orca.orca_haloed_global_grid_index(i,j));
            }
          }
        }
      }
      int total_is_node =
          std::count(local_orca.is_node.begin(), local_orca.is_node.end(), true);
      int total_is_ghost =
          std::count(local_orca.is_ghost.begin(), local_orca.is_ghost.end(), true);
      EXPECT(total_is_node <= indices.size());
      if (distributionName == "serial")
        EXPECT(local_orca.nb_used_nodes() == local_orca.nx() * local_orca.ny());
      EXPECT(total_is_ghost >= local_orca.nb_used_ghost_nodes());
      EXPECT(indices.size() == local_orca.nx() * local_orca.ny());

      {
        // diagnostics
        auto total_on_partition =
            std::count(this_partition.begin(), this_partition.end(), true);
        auto not_on_partition =
            std::count(this_partition.begin(), this_partition.end(), false);

        std::cout << "[" << cfg.mypart << "] grid.haloWest() " << grid.haloWest()
                  << " grid.haloEast() " << grid.haloEast()
                  << " grid.haloNorth() " << grid.haloNorth()
                  << " grid.haloSouth() " << grid.haloSouth()
                  << std::endl;
        std::cout << "[" << cfg.mypart << "]"
                  << " ix_orca_min " << local_orca.ix_min() << " ix_orca_max "
                  << local_orca.ix_max() << " iy_orca_min " << local_orca.iy_min()
                  << " iy_orca_max " << local_orca.iy_max() << " indices.size() "
                  << indices.size() << " nx * ny " << local_orca.nx() << " * " << local_orca.ny()
                  << " " << local_orca.nx() * local_orca.ny()
                  << " number on this partition " << total_on_partition
                  << " number not on partition " << not_on_partition << std::endl;

        output::Gmsh gmsh(std::string("surroundingRect") +
                              std::to_string(cfg.nparts) + "_" + gridname + "_" +
                              distributionName + "_" + std::to_string(halo) +
                              ".msh",
                          Config("coordinates", "xy") | Config("info", true));
      }

      const idx_t cell_width = 1;
      // Values below are for ORCA2_T grid.
      if (cfg.nparts == 2) {
        if (cfg.mypart == 0) {
          EXPECT(local_orca.iy_min() == -grid.haloSouth());
          EXPECT(local_orca.iy_max() == grid.ny() + halo);
          EXPECT(local_orca.ix_min() == -std::max(grid.haloWest(), halo));
          EXPECT(local_orca.ix_max() == 90 + halo);
          EXPECT(local_orca.nx() == 91 + halo + std::max(grid.haloWest(), halo));
          EXPECT(local_orca.ny() == grid.ny() + grid.haloSouth() + halo + 1);
        }
        if (cfg.mypart == 1) {
          EXPECT(local_orca.iy_min() == -grid.haloSouth());
          EXPECT(local_orca.iy_max() == grid.ny() + halo);
          EXPECT(local_orca.ix_min() == 90 - halo);
          EXPECT(local_orca.ix_max() == grid.nx() + halo);
          EXPECT(local_orca.nx() == grid.nx() + 2 * halo - 89);
          EXPECT(local_orca.ny() == grid.ny() + grid.haloSouth() + halo + 1);
        }
      }
    }
  }
}

} // namespace test
} // namespace atlas

int main(int argc, char **argv) { return atlas::test::run(argc, argv); }
