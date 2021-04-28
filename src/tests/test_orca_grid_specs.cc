/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/SpecRegistry.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test spec" ) {
    auto registered = [&]( const std::string& uid_or_name ) { return grid::SpecRegistry::has( uid_or_name ); };
    auto orca_name  = [&]( const std::string& uid_or_name ) {
        return grid::SpecRegistry::get( uid_or_name ).getString( "orca_name" );
    };
    auto orca_arrangement = [&]( const std::string& uid_or_name ) {
        return grid::SpecRegistry::get( uid_or_name ).getString( "orca_arrangement" );
    };
    auto dimensions = [&]( const std::string& uid_or_name ) {
        std::vector<int> v;
        grid::SpecRegistry::get( uid_or_name ).get( "dimensions", v );
        return v;
    };
    auto uid = [&]( const std::string& uid_or_name ) {
        return grid::SpecRegistry::get( uid_or_name ).getString( "uid" );
    };


    std::vector<std::string> uids;
    std::map<std::string, std::string> check_name;

    std::vector<std::string> grids{"ORCA2", "ORCA1", "eORCA1", "ORCA025", "eORCA025", "ORCA12", "eORCA12"};
    for ( const auto& grid : grids ) {
        std::vector<std::string> arrangements{"F", "T", "U", "V", "W"};
        for ( const auto& P : arrangements ) {
            std::string name = grid + "_" + P;
            EXPECT( registered( name ) );
            EXPECT( orca_name( name ) == grid );
            EXPECT( orca_arrangement( name ) == P );
            EXPECT_NO_THROW( uid( name ) );
            EXPECT_NO_THROW( dimensions( name ) );
            Log::info() << std::setw( 11 ) << std::left << name << "    " << uid( name ) << "    " << dimensions( name )
                        << std::endl;
            uids.push_back( uid( name ) );
            check_name[uids.back()] = name;
        }
    }

    for ( const auto& grid : uids ) {
        EXPECT( registered( grid ) );
        EXPECT( grid == uid( grid ) );
        EXPECT( orca_name( grid ) + "_" + orca_arrangement( grid ) == check_name[grid] );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
