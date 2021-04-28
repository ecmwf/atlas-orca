/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas-orca/Library.h"

#include "eckit/filesystem/PathName.h"

#include <string>
#include <vector>

namespace atlas {
namespace orca {

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
