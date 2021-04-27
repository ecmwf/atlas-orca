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

#include <cstdlib>

#include "eckit/filesystem/PathName.h"

#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Exception.h"
#include "atlas/io/FileStream.h"


#include "eckit/io/URLHandle.h"
#include "eckit/log/Bytes.h"


namespace atlas {
namespace orca {

struct AutoIndent {
    AutoIndent() {
        Log::info().indent();
    }
    ~AutoIndent() {
        if( Log::info() ) {
            Log::info().unindent();
        }
    }
};


inline eckit::Length curl_download( const std::string& url, const eckit::PathName& path ) {

    // This version manages to overcome redirects, which URLHandle can currently not cope with
    std::string command = "curl --progress-bar -L " + url + " --output " + path;
    Log::debug() << "+ " << command << std::endl;
    if( std::system( command.c_str() ) == 0 ) {
        return path.size();
    }
    return 0;
}


inline size_t download( const std::string& url, const eckit::PathName& path ) {
    atlas::Trace trace(Here(),"Downloading ORCA grid data");
    auto parent_directory = path.dirName();
    parent_directory.mkdir();
    ATLAS_ASSERT( parent_directory.exists() );
    Log::info() << "Downloading " << url << " to " << path << " ..." << std::endl;
    AutoIndent indent;
    eckit::PathName path_tmp = path+".download";
    auto length = eckit::URLHandle(url).saveInto(path_tmp);
    if( length <= 0 ) {
        path_tmp.unlink(true);
        return 0;
    }
    if( length < eckit::Length(10*1024) ) {
        io::InputFileStream file(path_tmp);
        std::string content;
        content.resize( path_tmp.size() );
        file.read( const_cast<char*>(content.data()),content.size());
        if( content.find("Error 404") ) {
            path_tmp.unlink(true);
            return 0;
        }
    }
    eckit::PathName::rename( path_tmp, path);
    trace.stop();
    Log::info() << "Download of " << eckit::Bytes(length) << " took " << trace.elapsed() << " s." << std::endl;
    return length;
};


}  // namespace orca
}  // namespace atlas
