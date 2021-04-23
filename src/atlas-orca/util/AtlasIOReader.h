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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>


#include "eckit/filesystem/PathName.h"
#include "eckit/log/ProgressTimer.h"

#include "atlas/runtime/Exception.h"
#include "atlas/io/atlas-io.h"

#include "atlas-orca/util/OrcaData.h"
#include "atlas-orca/util/OrcaDataFile.h"

namespace atlas {
namespace orca {

class AtlasIOReader {
public:
    AtlasIOReader( const util::Config& config = util::NoConfig() ) {
    }

    void read( const std::string& uri, OrcaData& data ) {

        OrcaDataFile file{uri};

        io::RecordReader reader(file);

        int version;
        reader.read("version",version).wait();
        if( version == 0 ) {
            reader.read("dimensions",data.dimensions);
            reader.read("pivot",data.pivot);
            reader.read("halo",data.halo);
            reader.read("longitude",data.lon);
            reader.read("latitude",data.lat);
            reader.read("flags",data.flags);
            reader.wait();
        }
        else {
            ATLAS_THROW_EXCEPTION("Unsupported version " << version );
        }
    }
};

} // namespace orca
}  // namespace atlas

//------------------------------------------------------------------------------------------------------
