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

#include "atlas-orca/Library.h"

#include <sstream>
#include <vector>

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/URI.h"
#include "eckit/utils/Tokenizer.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "atlas-orca/util/ComputeCachedPath.h"
#include "atlas-orca/util/Download.h"


namespace atlas {
namespace orca {

static std::vector<std::string> known_urls() {
    std::vector<std::string> urls;
    urls.emplace_back("http://get.ecmwf.int/atlas/grids/orca");
    urls.emplace_back("https://get.ecmwf.int/atlas/grids/orca");
    return urls;
}

static std::vector<std::string> search_paths() {
    std::vector<std::string> paths;
    eckit::Tokenizer tokenize(":");
    std::vector<std::string> tokenized;
    tokenize( Library::instance().dataPath(), tokenized);
    for( auto& t : tokenized ) {
        if( not t.empty() ) {
            paths.push_back(t);
        }
    }
    if( Library::instance().caching() ) {
        paths.push_back( Library::instance().cachePath() );
    }
    return paths;
}


class OrcaDataFile {
public:
    OrcaDataFile( const std::string& uri, const std::string& checksum = "" ) {
        uri_      = eckit::URI( uri );
        checksum_ = checksum;
        auto file_search = FileSearch{search_paths(),known_urls()};
        if ( uri_.scheme().find( "http" ) == 0 ) {
            std::string url = uri_.asRawString();

            Log::debug() << "Looking for " << file_search.file(url) << " in " << file_search.searchPath() << std::endl;
            int found = 0;
            if ( mpi::comm().rank() == 0 ) {
                found = file_search( url, path_ );
            }
            mpi::comm().broadcast(found,0);
            if( not found ) {
                Log::debug() << "File " << uri << " has not been found in " << file_search.searchPath() << std::endl;
            }
            if( found ) {
                if( mpi::comm().rank() != 0 ) {
                    // Find path also on non-zero ranks
                    ATLAS_ASSERT( file_search( url, path_ ) );
                }
                Log::debug() << "File " << uri << " has been found: " << path_ << std::endl;
            }
            else if ( Library::instance().caching() ) {
                path_ = ComputeCachedPath{known_urls()}(url);

                Log::debug() << "Caching enabled. Downloading " << uri << " to " << path_ << std::endl;

                if ( mpi::comm().rank() == 0 ) {
                   if ( download( url, path_ ) == 0 ) {
                        std::stringstream errmsg;
                        errmsg << "Could not download file from url " << url;
                        ATLAS_THROW_EXCEPTION( errmsg.str() );
                    }
                    if ( not path_.exists() ) {
                        ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << path_ );
                    }
                }
                mpi::comm().barrier();
            }
            else {
                ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << file_search.file(url) << " in " << file_search.searchPath() );
            }
        }
        else {
            bool found;
            path_ = uri_.path();
            found = path_.exists();
            if( not found ) {
                found = file_search( uri_.path(), path_ );
                if( found ) {
                    Log::debug() << "File " << uri << " has been found: " << path_ << std::endl;
                }
            }
            if( not found ) {
                ATLAS_THROW_EXCEPTION( "Could not locate orca grid data file " << uri );
            }
        }
    }

    const eckit::PathName& path() const { return path_; }

    operator const eckit::PathName&() const { return path(); }
    operator std::string() const { return path_.asString(); }
    operator const char*() const { return c_str(); }
    const char* c_str() const { return path_.localPath(); }

private:
    eckit::URI uri_;
    eckit::PathName path_;
    std::string checksum_;
};

}  // namespace orca
}  // namespace atlas
