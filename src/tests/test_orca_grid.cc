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
#include <fstream>
#include <iomanip>

#include "eckit/log/Bytes.h"

#include "atlas/grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/util/Config.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"
#include "atlas-orca/util/OrcaPeriodicity.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;


namespace atlas::test {

//-----------------------------------------------------------------------------

CASE( "test orca grid iterator" ) {
    struct Section {
        std::string gridname;
        size_t size;
    };

    std::vector<Section> sections{
        { "ORCA2_T", 27118 }, { "eORCA1_T", 120184 },
        //{"ORCA025_T", 1472282},
    };
    for ( auto& section : sections ) {
        std::string gridname = section.gridname;
        SECTION( gridname ) {
            OrcaGrid grid( gridname );

            EXPECT_EQ( grid.size(), section.size );

            Log::info() << "grid.footprint() = " << eckit::Bytes( static_cast<double>( grid.footprint() ) )
                        << std::endl;

            idx_t n = 0;
            {
                auto trace = Trace( Here(), "iterating" );
                for ( const auto& p : grid.lonlat() ) {
                    ++n;
                }
                trace.stop();
                Log::info() << "iterating took " << trace.elapsed() << " seconds" << std::endl;
            }
            EXPECT_EQ( n, grid.size() );
            Log::info() << "First point: " << grid.lonlat().front() << std::endl;
            Log::info() << "Last point: " << grid.lonlat().back() << std::endl;
        }
    }
}

CASE("test matchup between orca and regular ij indexing ") {
  std::vector<std::string> gridnames = {"ORCA2_T", "eORCA025_T"};

  for (const std::string& gridname : gridnames) {
    auto mypart = 0;
    auto orca_grid = OrcaGrid(gridname);

    StructuredGrid::YSpace yspace{grid::LinearSpacing{
        {-80., 90.}, orca_grid.ny(), true}};
    StructuredGrid::XSpace xspace{
        grid::LinearSpacing{{0., 360.}, orca_grid.nx(), false}};
    StructuredGrid regular_grid{xspace, yspace};

    SECTION("is the regular grid ij index unique?") {
      const idx_t size = regular_grid.size();
      std::set<gidx_t> ij_uid;
      for (gidx_t node = 0; node < size; ++node) {
        idx_t i, j;
        regular_grid.index2ij(node, i, j);
        ij_uid.insert(i*10*size + j);
      }

      if (size != ij_uid.size())
        std::cout << "[" << mypart
                  << "] number of duplicate regular grid ij UID points "
                  << size - ij_uid.size()
                  << "/" << size << std::endl;
      EXPECT(size == ij_uid.size());
    }

    SECTION("are the orca grid internal ij indices unique?") {
      const idx_t size = orca_grid.size();
      std::set<gidx_t> ij_uid;
      for (gidx_t node = 0; node < size; ++node) {
        idx_t i, j;
        orca_grid.index2ij(node, i, j);
        if (i >= orca_grid.nx() || i < 0) continue;
        if (j >= orca_grid.ny() || j < 0) continue;
        ij_uid.insert(i*10*size + j);
      }
      auto internal_size = orca_grid.nx()*orca_grid.ny();
      if (internal_size != ij_uid.size())
        std::cout << "[" <<  mypart
                  << "] number of duplicate orca grid ij UID points "
                  << internal_size - ij_uid.size()
                  << "/" << internal_size << std::endl;
      EXPECT(internal_size == ij_uid.size());
    }

    auto nx_orca_halo = orca_grid.haloWest() + orca_grid.nx() + orca_grid.haloEast(); 
    auto ny_orca_halo = orca_grid.haloSouth() + orca_grid.ny() + orca_grid.haloNorth();

    const std::array<std::int32_t, 2> dimensions{nx_orca_halo, ny_orca_halo};
    const std::array<std::int32_t, 4> halo{orca_grid.haloNorth(), orca_grid.haloWest(),
                                           orca_grid.haloSouth(), orca_grid.haloEast()};
    const std::array<double, 2> pivot{orca_grid.nx()/2 + 1, orca_grid.ny()};

    if (gridname == "ORCA2_T") {
      EXPECT(nx_orca_halo == 182);
      EXPECT(ny_orca_halo == 149);
      EXPECT(orca_grid.haloNorth() == 1);
      EXPECT(orca_grid.haloWest() == 1);
      EXPECT(orca_grid.haloSouth() == 1);
      EXPECT(orca_grid.haloEast() == 1);
      std::cout << "pivot[0]: " << pivot[0] << " =? " << 91 << std::endl;
      std::cout << "pivot[1]: " << pivot[1] << " =? " << 147 << std::endl;
      EXPECT(pivot[0] == 91);
      EXPECT(pivot[1] == 147);
    }

    atlas::orca::OrcaPeriodicity periodicity(dimensions, halo, pivot);

    SECTION("ORCA periodicity satisfies symmetries") {

      for (idx_t i = 0; i < dimensions[0]; ++i) {
        for (idx_t j = 0; j < dimensions[1]; ++j) {
          if (j == orca_grid.ny() && i > orca_grid.nx()/2) {
            auto master = periodicity(i, j);
            EXPECT(master.j == orca_grid.ny());
            EXPECT(master.i == orca_grid.nx() + 2 - i);
          }
          // top left corner there are two nodes that don't make sense to me
          if (j == orca_grid.ny() + 1 && i > 2) {
            auto master = periodicity(i, j);
            EXPECT(master.j == orca_grid.ny() - 1);
            EXPECT(master.i == orca_grid.nx() + 2 - i);
          }
        }
      }
    }

    SECTION("Are the boundary symmetries present in the orca grid ij indices?") {
      const idx_t size = orca_grid.size();
      atlas::PointLonLat lonlat, lonlat_halo;
      for (gidx_t node = 0; node < size; ++node) {
        idx_t i, j;
        orca_grid.index2ij(node, i, j);
        auto pivot = orca_grid.nx() / 2 + 1;
        if (gridname.back() == 'U') {
          pivot = orca_grid.nx() / 2;
        }
        if (j == orca_grid.ny() && i > 2) {
          // check northfold boundary
          //     - second row is swapped version of second from top row
          lonlat = orca_grid.lonlat(i, j);
          lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j - 2);
          EXPECT(orca_grid.periodicIndex(i, j) ==
                 orca_grid.periodicIndex(orca_grid.nx() - i, j - 2));
          EXPECT(lonlat[0] == lonlat_halo[0]);
          EXPECT(lonlat[1] == lonlat_halo[1]);
        }
        if (j == orca_grid.ny() - 1) {
          // check northfold boundary - centre fold row is mirrored about the central pivot
          if (i <= pivot) continue; // halo points are right of pivot on this row.
          // count from right-hand-side of grid.
          lonlat = orca_grid.lonlat(i, j);
          lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j);
          EXPECT(orca_grid.periodicIndex(i, j) ==
                 orca_grid.periodicIndex(orca_grid.nx() - i, j));
          EXPECT(lonlat[0] == lonlat_halo[0]);
          EXPECT(lonlat[1] == lonlat_halo[1]);
        }
        if (i > orca_grid.nx()) {
          // check east-west boundary
          lonlat = orca_grid.lonlat(i, j);
          lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j);
          EXPECT(orca_grid.periodicIndex(i, j) ==
                 orca_grid.periodicIndex(orca_grid.nx() - i, j));
          EXPECT(lonlat[0] == lonlat_halo[0]);
          EXPECT(lonlat[1] == lonlat_halo[1]);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
