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


#include "eckit/filesystem/PathName.h"

#include "atlas/library.h"
#include "atlas/runtime/Exception.h"

#include "atlas-orca/Library.h"

#include <string>
#include <vector>

namespace atlas {
namespace orca {

class FileSearch {
public:
    FileSearch( const std::vector<std::string> search_paths, const std::vector<std::string>& known_urls ) :
        search_paths_{search_paths},
        known_urls_{known_urls} {
    }
    bool operator()( const std::string& url, eckit::PathName& path ) const {
        std::string search_file = file(url);
        for( auto& search_path : search_paths_ ) {
            std::string join = "";
            if( not search_path.empty() && search_path[search_path.size()-1] != '/' ) {
                join = "/";
            }
            eckit::PathName search = search_path + join + search_file;
            if( search.exists() ) {
                path = search;
                return true;
            }
        }
        return false;
    }

    std::string file( const std::string& url ) const {
        if ( url.rfind( "http", 0 ) == 0 ) {
            bool is_known_url = false;
            for ( const auto& known_url : known_urls_ ) {
                if ( url.find( known_url ) == 0 ) {
                    is_known_url         = true;
                    std::string filepath = url.substr( known_url.size() );
                    return "atlas/grids/orca" + filepath;
                }
            }
            if ( not is_known_url ) {
                std::string filename = url.substr( url.find_last_of( "/" ), url.size() );
                return "atlas/grids/orca/unknown" + filename;
            }
            ATLAS_THROW_EXCEPTION("Should not be here");
        }
        else {
            return url;
        }
    }

    std::string searchPath() const {
        std::stringstream joined;
        for( size_t i=0; i<search_paths_.size(); ++i ) {
            if( i > 0 ) {
                joined << ":";
            }
            joined << search_paths_[i];
        }
        return joined.str();
    }

private:
    std::vector<std::string> known_urls_;
    std::vector<std::string> search_paths_;
};

class ComputeCachedPath {
public:
    ComputeCachedPath( const std::vector<std::string>& known_urls ) : known_urls_{known_urls} {}
    eckit::PathName operator()( const std::string& url ) const {
        eckit::PathName path;
        bool is_known_url = false;
        for ( const auto& known_url : known_urls_ ) {
            if ( url.find( known_url ) == 0 ) {
                is_known_url         = true;
                std::string filepath = url.substr( known_url.size() );
                path                 = Library::instance().cachePath() + "/atlas/grids/orca" + filepath;
            }
        }
        if ( not is_known_url ) {
            std::string filename = url.substr( url.find_last_of( "/" ), url.size() );
            path                 = Library::instance().cachePath() + "/atlas/grids/orca/unknown" + filename;
        }
        return path;
    }

private:
    std::vector<std::string> known_urls_;
};

}  // namespace orca
}  // namespace atlas
