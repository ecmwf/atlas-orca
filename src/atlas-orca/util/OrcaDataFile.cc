/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas-orca/util/OrcaDataFile.h"

#include <sstream>
#include <vector>

#include "atlas/runtime/Log.h"
#include "eckit/utils/Tokenizer.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

#include "atlas-orca/Library.h"
#include "atlas-orca/util/ComputeCachedPath.h"
#include "atlas-orca/util/Download.h"


namespace atlas::orca {

OrcaDataFile::OrcaDataFile( const std::string& uri, const std::string& checksum ) : uri_( uri ), checksum_( checksum ) {
    static const std::vector<std::string> KNOWN_URLS{
        "http://get.ecmwf.int/repository/atlas/grids/orca",
        "https://get.ecmwf.int/repository/atlas/grids/orca",
    };

    static const std::vector<std::string> SEARCH_PATHS = []() {
        std::vector<std::string> paths;
        eckit::Tokenizer( ":", false )( Library::instance().dataPath(), paths );
        if ( Library::instance().caching() ) {
            paths.emplace_back( Library::instance().cachePath() );
        }
        return paths;
    }();

    auto file_search = FileSearch{ SEARCH_PATHS, KNOWN_URLS };

    if ( uri_.scheme().find( "http" ) == 0 ) {
        auto url = uri_.asRawString();

        Log::debug() << "Looking for " << file_search.file( url ) << " in " << file_search.searchPath() << std::endl;
        int found = 0;
        if ( mpi::comm().rank() == 0 ) {
            found = file_search( url, path_ ) ? 1 : 0;
        }
        mpi::comm().broadcast( found, 0 );
        if ( found == 0 ) {
            Log::debug() << "File " << uri << " has not been found in " << file_search.searchPath() << std::endl;
        }
        if ( found == 1 ) {
            if ( mpi::comm().rank() != 0 ) {
                // Find path also on non-zero ranks
                ATLAS_ASSERT( file_search( url, path_ ) );
            }
            Log::debug() << "File " << uri << " has been found: " << path_ << std::endl;
        }
        else if ( Library::instance().caching() ) {
            path_ = ComputeCachedPath{ KNOWN_URLS }( url );

            Log::debug() << "Caching enabled. Downloading " << uri << " to " << path_ << std::endl;

            if ( mpi::comm().rank() == 0 ) {
                if ( download( url, path_ ) == 0 ) {
                    ATLAS_THROW_EXCEPTION( "Could not download file from url " << url );
                }
                if ( not path_.exists() ) {
                    ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << path_ );
                }
            }
            mpi::comm().barrier();
        }
        else {
            ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << file_search.file( url ) << " in "
                                                                           << file_search.searchPath() );
        }
    }
    else {
        path_      = uri_.path();
        auto found = path_.exists();
        if ( not found ) {
            found = file_search( uri_.path(), path_ );
            if ( found ) {
                Log::debug() << "File " << uri << " has been found: " << path_ << std::endl;
            }
        }
        if ( not found ) {
            ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << uri );
        }
    }
}

}  // namespace atlas::orca
