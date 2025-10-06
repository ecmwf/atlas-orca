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
  std::vector<std::string> gridnames = {
    "ORCA2_T",
    // "eORCA025_T",
  };

  for (const std::string& gridname : gridnames) {
    SECTION(gridname + ": is the regular grid ij index unique?") {
      auto mypart = 0;
      auto orca_grid = OrcaGrid(gridname);

      StructuredGrid::YSpace yspace{grid::LinearSpacing{
          {-80., 90.}, orca_grid.ny(), true}};
      StructuredGrid::XSpace xspace{
          grid::LinearSpacing{{0., 360.}, orca_grid.nx(), false}};
      StructuredGrid regular_grid{xspace, yspace};

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

    SECTION(gridname + ": are the orca grid internal ij indices unique?") {
      auto mypart = 0;
      auto orca_grid = OrcaGrid(gridname);

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


    SECTION(gridname + ": ORCA periodicity satisfies symmetries") {
      auto orca_grid = OrcaGrid(gridname);

      auto nx_orca_halo = orca_grid.haloWest() + orca_grid.nx() + orca_grid.haloEast();
      auto ny_orca_halo = orca_grid.haloSouth() + orca_grid.ny() + orca_grid.haloNorth();

      const std::array<std::int32_t, 2> dimensions{static_cast<std::int32_t>(nx_orca_halo), static_cast<std::int32_t>(ny_orca_halo)};
      const std::array<std::int32_t, 4> halo{orca_grid.haloNorth(), orca_grid.haloWest(),
                                            orca_grid.haloSouth(), orca_grid.haloEast()};
      const std::array<double, 2> pivot{static_cast<double>(orca_grid.nx()/2 + 1), static_cast<double>(orca_grid.ny())};

      if (gridname == "ORCA2_T") {
        EXPECT_EQ(nx_orca_halo, 182);
        EXPECT_EQ(ny_orca_halo, 149);
        EXPECT_EQ(orca_grid.haloNorth(), 1);
        EXPECT_EQ(orca_grid.haloWest(), 1);
        EXPECT_EQ(orca_grid.haloSouth(), 1);
        EXPECT_EQ(orca_grid.haloEast(), 1);
        EXPECT_EQ(pivot[0], 91);
        EXPECT_EQ(pivot[1], 147);
      }
  
      atlas::orca::OrcaPeriodicity periodicity(dimensions, halo, pivot);

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

    SECTION(gridname + ": Are the boundary symmetries present in the orca grid ij indices?") {
      auto orca_grid = OrcaGrid(gridname);
      auto pivot = orca_grid.nx() / 2 + 1;
      if (gridname.back() == 'U') {
        pivot = orca_grid.nx() / 2;
      }

      const idx_t size = orca_grid.size();
      atlas::PointLonLat lonlat, lonlat_halo;
      for (gidx_t node = 0; node < size; ++node) {
        idx_t i, j;
        orca_grid.index2ij(node, i, j);
        if (gridname == "ORCA2_T") {
          if (j == orca_grid.ny() && i > 2) {
            // check northfold boundary
            //     - second row is swapped version of second from top row
            lonlat = orca_grid.lonlat(i, j);
            lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j - 2);
            EXPECT_EQ(orca_grid.periodicIndex(i, j),
                      orca_grid.periodicIndex(orca_grid.nx() - i, j - 2));
            EXPECT_EQ(lonlat[0], lonlat_halo[0]);
            EXPECT_EQ(lonlat[1], lonlat_halo[1]);
          }
          if (j == orca_grid.ny() - 1) {
            // check northfold boundary - centre fold row is mirrored about the central pivot
            if (i <= pivot) continue; // halo points are right of pivot on this row.
            // count from right-hand-side of grid.
            lonlat = orca_grid.lonlat(i, j);
            lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j);
            EXPECT_EQ(orca_grid.periodicIndex(i, j),
                      orca_grid.periodicIndex(orca_grid.nx() - i, j));
            EXPECT_EQ(lonlat[0], lonlat_halo[0]);
            EXPECT_EQ(lonlat[1], lonlat_halo[1]);
          }
        }
        if (i > orca_grid.nx()) {
          // check east-west boundary
          lonlat = orca_grid.lonlat(i, j);
          lonlat_halo = orca_grid.lonlat(orca_grid.nx() - i, j);
          EXPECT_EQ(orca_grid.periodicIndex(i, j),
                    orca_grid.periodicIndex(orca_grid.nx() - i, j));
          EXPECT_EQ(lonlat[0], lonlat_halo[0]);
          EXPECT_EQ(lonlat[1], lonlat_halo[1]);
        }
        if (i < 0) {
          // check east-west boundary
          lonlat_halo = orca_grid.lonlat(i, j);
          lonlat = orca_grid.lonlat(orca_grid.nx() + i, j);
          EXPECT_EQ(orca_grid.periodicIndex(i, j),
                    orca_grid.periodicIndex(orca_grid.nx() + i, j));
          EXPECT_EQ(lonlat[0], lonlat_halo[0]);
          EXPECT_EQ(lonlat[1], lonlat_halo[1]);
        }

      }
    }
  }
}

CASE("periodicity") {
  SECTION("ORCA2_T") {
    OrcaGrid g("ORCA2_T");
    EXPECT_EQ(g.nx(), 180);
    EXPECT_EQ(g.ny(), 147);
    EXPECT_EQ(g.periodicIndex(-1,-1), 180);
    EXPECT_EQ(g.periodicIndex(g.nx()-1,-1), 180);
    EXPECT_EQ(g.periodicIndex(0,-1), 1);
    EXPECT_EQ(g.periodicIndex(0,0), 183);
    EXPECT_EQ(g.periodicIndex(g.nx(),0), 183);
    EXPECT_EQ(g.periodicIndex(g.nx()-1,0), 362);
    EXPECT_EQ(g.periodicIndex(-1,0), 362);

    // now go outside of built-in halo, which should wrap around into core region points
    EXPECT_EQ(g.periodicIndex(-10,0), g.periodicIndex(-1,0) - 9);
    EXPECT_EQ(g.periodicIndex(-10,0), 353); // and evaluated
    EXPECT_EQ(g.periodicIndex(g.nx()+10,0), g.periodicIndex(0,0) + 10);
    EXPECT_EQ(g.periodicIndex(g.nx()+10,0), 193); // and evaluated

    // north fold for ORCA2_T already starts in row ny-1
    int ny = g.ny();
    EXPECT_EQ(g.periodicIndex(0,  ny-1), 26755); // Beginning of row
    EXPECT_EQ(g.periodicIndex(90, ny-1), 26845); // pivot center point
    EXPECT_EQ(g.periodicIndex(91, ny-1), 26844); // Just right of pivot, folding back
    EXPECT_EQ(g.periodicIndex(92, ny-1), 26843); // ..
    EXPECT_EQ(g.periodicIndex(179,ny-1), 26756); // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIndex(g.nx(),ny-1), 26755); // Periodic point of row must return beginning of row
    EXPECT_EQ(g.periodicIndex(-1,  ny-1), 26756);   // Periodic point must return last of row (179,ny-1), which folds back to second point of row

    // North halo
    EXPECT_EQ(g.periodicIJ(180, ny), g.periodicIJ(0,  ny));    // Periodic within fold
    EXPECT_EQ(g.periodicIJ(0,  ny), g.periodicIJ(180, ny-2));  // This also folds back but into halo in core rows
    EXPECT_EQ(g.periodicIJ(180, ny-2), g.periodicIJ(0, ny-2)); // And this is periodic with core point
    EXPECT_EQ(g.periodicIndex(0, ny-2), 26573);
    EXPECT_EQ(g.periodicIJ(-1,  ny), g.periodicIJ(179, ny));
    EXPECT_EQ(g.periodicIJ(179, ny), g.periodicIJ(1,ny-2));
    EXPECT_EQ(g.periodicIndex(1,ny-2), 26574);
    // extra halo:
    EXPECT_EQ(g.periodicIJ(180+10, ny), g.periodicIJ(0+10,  ny));       // Periodic within fold
    EXPECT_EQ(g.periodicIJ(0+10, ny),   g.periodicIJ(180-10,  ny-2));   // Folds back into core region
    EXPECT_EQ(g.periodicIndex(180+10, ny-2),26583);   

    // Extra North halo not part of grid definition
    EXPECT_EQ(g.periodicIJ(0, (ny-1)+10), g.periodicIJ(180,(ny-1)-10));
    EXPECT_EQ(g.periodicIJ(180, (ny-1)-10), g.periodicIJ(0,(ny-1)-10));
    EXPECT_EQ(g.periodicIndex(0,(ny-1)-10), 24935);
    EXPECT_EQ(g.periodicIndex(0,(ny-1)+10), 24935);
  }

  SECTION("ORCA1_T") {
    OrcaGrid g("ORCA1_T");
    EXPECT_EQ(g.nx(), 360);
    EXPECT_EQ(g.ny(), 290);
    EXPECT_EQ(g.periodicIndex(-1,-1), 360);
    EXPECT_EQ(g.periodicIndex(g.nx()-1,-1), 360);
    EXPECT_EQ(g.periodicIndex(0,-1), 1);
    EXPECT_EQ(g.periodicIndex(0,0), 363);
    EXPECT_EQ(g.periodicIndex(g.nx(),0), 363);
    EXPECT_EQ(g.periodicIndex(g.nx()-1,0), 722);
    EXPECT_EQ(g.periodicIndex(-1,0), 722);

    // now go outside of built-in halo, which should wrap around into core region points
    EXPECT_EQ(g.periodicIndex(-10,0), g.periodicIndex(-1,0) - 9);
    EXPECT_EQ(g.periodicIndex(-10,0), 713); // and evaluated
    EXPECT_EQ(g.periodicIndex(g.nx()+10,0), g.periodicIndex(0,0) + 10);
    EXPECT_EQ(g.periodicIndex(g.nx()+10,0), 373); // and evaluated

    // north fold for ORCA1_T does not start in row ny-1, but half increment above
    int ny = g.ny();
    int nx = g.nx();
    EXPECT_EQ(g.periodicIndex(0,   ny-1), 104981); // Beginning of row
    EXPECT_EQ(g.periodicIndex(180, ny-1), 105161); // pivot center point
    EXPECT_EQ(g.periodicIndex(181, ny-1), 105162); // Just right of pivot, folding back
    EXPECT_EQ(g.periodicIndex(182, ny-1), 105163); // ..
    EXPECT_EQ(g.periodicIndex(359, ny-1), 105340);  // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIndex(g.nx(),ny-1), 104981); // Periodic point of row must return beginning of row
    EXPECT_EQ(g.periodicIndex(-1,  ny-1), 105340);   // Periodic point must return last of row (179,ny-1), which folds back to second point of row

    // North halo folds
    EXPECT_EQ(g.periodicIndex(0,   ny), 105340); // Beginning of row folds to last point of owned (nx-1,ny-1)
    EXPECT_EQ(g.periodicIndex(1,   ny), 105339); // Second of row folds to (nx-1 -1,ny-1)
    EXPECT_EQ(g.periodicIndex(2,   ny), 105338); // Third of row folds to (nx-1 -2,ny-1)
    EXPECT_EQ(g.periodicIndex(180, ny), 105160); // Just left of pivot center point
    EXPECT_EQ(g.periodicIndex(181, ny), 105159); // Just right of pivot center point
    EXPECT_EQ(g.periodicIndex(182, ny), 105158); // ..
    EXPECT_EQ(g.periodicIndex(359, ny), 104981);  // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIJ(359, ny), g.periodicIJ(0, ny-1));  // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIJ(nx, ny), g.periodicIJ(0, ny));  // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIndex(nx,ny), 105340); // Periodic point of row must return beginning of row
    EXPECT_EQ(g.periodicIJ(-1, ny), g.periodicIJ(nx-1, ny));  // Last of row folds back to second point of row
    EXPECT_EQ(g.periodicIndex(-1,  ny), 104981);   // Periodic point must return last of row (179,ny-1), which folds back to second point of row
  }


}

//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
