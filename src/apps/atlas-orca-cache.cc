/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <iostream>
#include <set>
#include <string>
#include <vector>

#if ATLAS_ORCA_HAVE_ECKIT_CODEC
#include "eckit/codec/codec.h"
#else
// Backward compatibility, DEPRECATED!
#include "atlas/io/atlas-io.h"
namespace eckit::codec {
using RecordReader  = atlas::io::RecordReader;
using InvalidRecord = atlas::io::InvalidRecord;
}  // namespace eckit::codec
#endif

#include "atlas/runtime/AtlasTool.h"

#include "atlas-orca/Library.h"
#include "atlas-orca/util/ComputeCachedPath.h"
#include "atlas-orca/util/Download.h"


namespace atlas::orca {

//----------------------------------------------------------------------------------------------------------------------

struct Tool : public atlas::AtlasTool {
    bool serial() override { return true; }
    int execute( const Args& ) override;
    std::string briefDescription() override { return "Create binary grid data files "; }
    std::string usage() override { return name() + " <file> --grid=NAME [OPTION]... [--help,-h]"; }
    std::string longDescription() override {
        return "Create binary grid data files \n"
               "\n"
               "       <file>: input file";
    }

    Tool( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "grid", "Grid to cache" ) );
        add_option( new SimpleOption<bool>( "overwrite", "Overwrite existing cache" ) );
    }
};

//------------------------------------------------------------------------------------------------------

int Tool::execute( const Args& args ) {
    auto specs = util::Config( Library::instance().gridsPath() );

    std::string grid = "all";
    args.get( "grid", grid );

    std::vector<std::string> grids;
    if ( grid != "all" ) {
        if ( not specs.has( grid ) ) {
            Log::error() << "Grid " << grid << " is not a known atlas-orca grid" << std::endl;
            return failed();
        }
        grids.push_back( grid );
    }
    else {
        for ( auto& id : specs.keys() ) {
            grids.push_back( id );
        }
    }
    std::set<std::string> urls;
    for ( auto& id : grids ) {
        auto spec = specs.getSubConfiguration( id );
        if ( spec.has( "data" ) ) {
            std::string uri = spec.getString( "data" );
            if ( uri.find( "http" ) == 0 ) {
                urls.insert( uri );
            }
        }
    }
    ComputeCachedPath compute_cached_path(
        { "https://get.ecmwf.int/repository/atlas/grids/orca", "http://get.ecmwf.int/repository/atlas/grids/orca" } );

    std::vector<std::string> failed_urls;
    for ( const auto& url : urls ) {
        auto path = compute_cached_path( url );
        if ( path.exists() ) {
            Log::info() << "File " << url << " was already found in cache: " << path << std::endl;
        }
        if ( !path.exists() || args.getBool( "overwrite", false ) ) {
            if ( download( url, path ) == 0 ) {
                Log::error() << "Error: Could not download file from url " << url << "\n" << std::endl;
                failed_urls.push_back( url );
            }
            else {
                try {
                    eckit::codec::RecordReader record( path );
                    auto metadata = record.metadata( "dimensions" );
                }
                catch ( const eckit::codec::InvalidRecord& ) {
                    Log::error() << "Error: Downloaded file " << path << " is invalid. Deleting...\n" << std::endl;
                    path.unlink( true );
                    failed_urls.push_back( url );
                }
            }
        }
    }
    if ( not failed_urls.empty() ) {
        Log::error() << "\nErrors occured while downloading from following urls:" << std::endl;
        for ( auto& url : failed_urls ) {
            Log::error() << "    " << url << std::endl;
        }
        return failed();
    }
    return success();
}

}  // namespace atlas::orca


//------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    atlas::orca::Tool tool( argc, argv );
    return tool.start();
}
