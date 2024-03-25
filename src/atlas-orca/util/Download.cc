/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Download.h"

#include "atlas-orca/Library.h"

#include <cstdlib>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/io/FileHandle.h"
#include "eckit/log/Bytes.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#ifdef ATLAS_ORCA_HAVE_ECKIT_URLHANDLE
#include "eckit/io/URLHandle.h"
#endif


namespace atlas::orca {

struct AutoIndent {
    AutoIndent() { Log::info().indent(); }

    AutoIndent( const AutoIndent& ) = delete;
    AutoIndent( AutoIndent&& )      = delete;

    void operator=( const AutoIndent& ) = delete;
    void operator=( AutoIndent&& )      = delete;

    ~AutoIndent() {
        if ( Log::info() ) {
            Log::info().unindent();
        }
    }
};


inline eckit::Length curl_download( const std::string& url, const eckit::PathName& path ) {
    // Fallback for when URLHandle cannot handle
    std::stringstream command;
    command << "curl -sS -k -L " << url << " --output " << path;
    Log::debug() << "+ " << command.str() << std::endl;
    auto errcode = std::system( command.str().c_str() );
    if ( errcode == 0 ) {
        return path.size();
    }
    return 0;
}


size_t download( const std::string& url, const eckit::PathName& path ) {
    atlas::Trace trace( Here(), "Downloading ORCA grid data" );
    auto parent_directory = path.dirName();
    parent_directory.mkdir();
    ATLAS_ASSERT( parent_directory.exists() );
    Log::info() << "Downloading " << url << " to " << path << " ..." << std::endl;
    AutoIndent indent;
    eckit::PathName path_tmp = path + ".download";
    eckit::Length length;
#ifdef ATLAS_ORCA_HAVE_ECKIT_URLHANDLE
    try {
        length = eckit::URLHandle( url ).saveInto( path_tmp );
        if ( length <= 0 && eckit_version_int() <= 11601 /*1.16.1*/ ) {
            // Problems with eckit::URLHandle fixed in further version
            Log::warning() << "Download failed with eckit::URLHandle. Trying again with curl system call." << std::endl;
            length = curl_download( url, path_tmp );
        }
    }
    catch ( eckit::SeriousBug ) {
        Log::warning() << "Download failed with eckit::URLHandle. Trying again with curl system call." << std::endl;
#else
    try {
#endif
        length = curl_download( url, path_tmp );
    }
    catch ( ... ) {
        length = 0;
    }

    if ( length <= 0 ) {
        if ( path_tmp.exists() ) {
            path_tmp.unlink( true );
        }
        return 0;
    }
    if ( length < eckit::Length( 10 * 1024 ) ) {
        std::string content;
        content.resize( path_tmp.size() );

        eckit::FileHandle file( path_tmp );
        file.openForRead();
        file.read( const_cast<char*>( content.data() ), static_cast<long>( content.size() ) );
        file.close();

        if ( content.find( "Error 404" ) != 0 ) {
            path_tmp.unlink( true );
            return 0;
        }
    }
    eckit::PathName::rename( path_tmp, path );
    trace.stop();
    Log::info() << "Download of " << eckit::Bytes( static_cast<double>( length ) ) << " took " << trace.elapsed()
                << " s." << std::endl;
    return length;
};


}  // namespace atlas::orca
