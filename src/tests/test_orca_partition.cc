/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/system/ResourceUsage.h"

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "atlas/util/Geometry.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "atlas-orca/util/PointIJ.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test orca partition is full" ) {
    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        "eORCA1_T",  //
        //"eORCA025_T",  //
    };
    auto methodnames = std::vector<std::string>{
        "bands",   //
        "equal_regions",  //
        "checkerboard",  //
    };
    for ( auto gridname : gridnames ) {
      auto grid = Grid( gridname );
      for ( auto method : methodnames ) {
        SECTION( gridname + std::string(" ") + method) {
            idx_t mypart = static_cast<idx_t>(mpi::rank());
            auto meshgenerator = MeshGenerator{"orca"};
            grid::Partitioner partitioner(method, atlas::mpi::size());
            auto mesh = meshgenerator.generate( grid, partitioner );
            auto part = array::make_view<int, 1>(mesh.nodes().partition());

            for (int p = 0; p<part.size(); p++) {
              if (part(p) < 0)
                std::cout << "[" << mypart << "] part(" << p << ") " << part(p) << std::endl;
              ASSERT(part(p) >= 0);
            }
        }
      }
    }
}

CASE( "test orca partition construction inside meshgenerator.generate" ) {

    auto gridnames = std::vector<std::string>{
        "ORCA2_T",   //
        "eORCA1_T",  //
        //"eORCA025_T",  //
    };
    auto methodnames = std::vector<std::string>{
        "bands",   //
        "equal_regions",  //
        "checkerboard",  //
    };
    for ( auto gridname : gridnames ) {
      auto grid = Grid( gridname );
      for ( auto method : methodnames ) {
        SECTION( gridname + std::string(" ") + method) {
          auto meshgenerator = MeshGenerator{"orca"};
          grid::Partitioner partitioner(method, atlas::mpi::size());
          auto mesh = meshgenerator.generate(grid, partitioner);
        }
      }
    }
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
