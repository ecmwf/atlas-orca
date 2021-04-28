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

#include <sstream>
#include <vector>

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/URI.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "atlas-orca/util/ComputeCachedPath.h"
#include "atlas-orca/util/Download.h"


namespace atlas {
namespace orca {

class OrcaDataFile {
public:
    OrcaDataFile( const std::string& uri, const std::string& checksum = "" ) {
        uri_      = eckit::URI( uri );
        checksum_ = checksum;
        if ( uri_.scheme().find( "http" ) == 0 ) {
            std::vector<std::string> known_urls{
                uri_.scheme() + "://get.ecmwf.int/atlas/grids/orca",
            };
            std::string url = uri_.asRawString();
            path_           = ComputeCachedPath{known_urls}( url );

            if ( mpi::comm().rank() == 0 ) {
                if ( not path_.exists() ) {
                    if ( Library::instance().download() ) {
                        if ( download( url, path_ ) == 0 ) {
                            std::stringstream errmsg;
                            errmsg << "Could not download file from url " << url;
                            ATLAS_THROW_EXCEPTION( errmsg.str() );
                        }
                    }
                    else {
                        std::stringstream errmsg;
                        errmsg << "File " << uri << " has not been found in cache: " << path_;
                        errmsg << "\nPlease run 'atlas-orca-cache' to populate cache, or export environment variable "
                                  "ATLAS_ORCA_DOWNLOAD=1";
                        ATLAS_THROW_EXCEPTION( errmsg.str() );
                    }
                }
                else {
                    Log::debug() << "File " << uri << " has already been found in cache: " << path_ << std::endl;
                }
            }
        }
        else {
            path_ = uri_.path();
        }
        if ( not path_.exists() ) {
            ATLAS_THROW_EXCEPTION( "Could not locate orca grid file " << path_ );
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
