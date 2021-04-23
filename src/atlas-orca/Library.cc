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

#include "eckit/filesystem/PathName.h"
#include "eckit/config/Resource.h"

#include "atlas/library/Library.h"
#include "atlas/util/Spec.h"

#include "atlas-orca/Library.h"
#include "atlas-orca/version.h"


namespace atlas {
class Grid;
}


namespace atlas {
namespace orca {


REGISTER_LIBRARY( Library );


Library::Library() : Plugin( "atlas-orca" ) {}


const Library& Library::instance() {
    static Library library;
    return library;
}


std::string Library::version() const {
    return atlas_orca_version();
}


std::string Library::gitsha1( unsigned int count ) const {
    std::string sha1 = atlas_orca_git_sha1();
    return sha1.empty() ? "not available" : sha1.substr( 0, std::min( count, 40U ) );
}


void Library::init() {
    Plugin::init();

    auto specs = util::Spec( grids() );
    for ( auto& id : specs.keys() ) {
        util::SpecRegistry<Grid>::enregister( id, specs.getSubConfiguration( id ) );
    }
}

std::string Library::cachePath() const {
    static std::string tmpdir = eckit::Resource<std::string>("$TMPDIR","/tmp");
    static std::string ATLAS_CACHE_PATH =
            eckit::PathName(eckit::LibResource<eckit::PathName, atlas::Library>(
                                "atlas-cache-path;$ATLAS_CACHE_PATH",tmpdir+"/cache"));
    static std::string ATLAS_ORCA_CACHE_PATH =
            eckit::PathName(eckit::LibResource<eckit::PathName, atlas::orca::Library>(
                                "atlas-orca-cache-path;$ATLAS_ORCA_CACHE_PATH",ATLAS_CACHE_PATH));
    return ATLAS_ORCA_CACHE_PATH;
}

bool Library::download() const {
    static bool ATLAS_ORCA_CACHE_DOWNLOAD =
            bool(eckit::LibResource<bool, atlas::orca::Library>(
                                "atlas-orca-download;$ATLAS_ORCA_DOWNLOAD",false));
    return ATLAS_ORCA_CACHE_DOWNLOAD;
}

std::string Library::grids() const {
    static std::string ATLAS_ORCA_GRIDS = eckit::LibResource<std::string, atlas::orca::Library>(
                           "atlas-orca-grids;$ATLAS_ORCA_GRIDS","~atlas-orca/etc/atlas-orca/grids.yaml");
    return ATLAS_ORCA_GRIDS;
}

}  // namespace orca
}  // namespace atlas
