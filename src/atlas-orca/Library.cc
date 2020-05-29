/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <string>

#include "atlas-orca/Library.h"
#include "atlas-orca/atlas_orca_config.h"


namespace atlas {
namespace orca {


REGISTER_LIBRARY( Library );


Library::Library() : eckit::system::Library( "atlas-orca" ) {}


const Library& Library::instance() {
    static Library library;
    return library;
}


const void* Library::addr() const {
    return this;
}


std::string Library::version() const {
    return atlas_orca_version();
}


std::string Library::gitsha1( unsigned int count ) const {
    std::string sha1 = atlas_orca_git_sha1();
    return sha1.empty() ? "not available" : sha1.substr( 0, std::min( count, 40U ) );
}


}  // namespace orca
}  // namespace atlas
