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

#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "atlas-orca/grid/OrcaGrid.h"

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

//-----------------------------------------------------------------------------

}  // namespace atlas::test


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
