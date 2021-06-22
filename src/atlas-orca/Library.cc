/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "eckit/config/Resource.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/runtime/Main.h"
#include "atlas-orca/Library.h"
#include "atlas-orca/version.h"
#include "atlas/grid/SpecRegistry.h"
#include "atlas/library/Library.h"


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
    auto grids = util::Config( gridsPath() );
    for ( auto& id : grids.keys() ) {
        Log::debug() << "Plugin atlas-orca registering grid " << id << std::endl;
        grid::SpecRegistry::add( id, grids.getSubConfiguration( id ) );
    }
}

bool Library::caching() const {
    static bool ATLAS_ORCA_CACHING =
        bool( eckit::LibResource<bool, atlas::orca::Library>( "atlas-orca-caching;$ATLAS_ORCA_CACHING", false ) );
    return ATLAS_ORCA_CACHING;
}

std::string Library::dataPath() const {
    return atlas::Library::instance().dataPath()+":~atlas-orca/share";
}

std::string Library::cachePath() const {
    return atlas::Library::instance().cachePath();
}

std::string Library::gridsPath() const {
    static std::string ATLAS_ORCA_GRIDS_PATH = eckit::LibResource<std::string, atlas::orca::Library>(
        "atlas-orca-grids-path;$ATLAS_ORCA_GRIDS_PATH", "" );
    ATLAS_ASSERT( not ATLAS_ORCA_GRIDS_PATH.empty() );
    return ATLAS_ORCA_GRIDS_PATH;
}

}  // namespace orca
}  // namespace atlas
