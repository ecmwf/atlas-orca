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

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "eckit/filesystem/PathName.h"

#include "atlas/runtime/Exception.h"

#include "atlas-orca/util/OrcaData.h"
#include "atlas-orca/util/OrcaDataFile.h"


namespace atlas::orca {

class NetCDFReader {
public:
    explicit NetCDFReader( const util::Config& = util::NoConfig() );

    void read( const std::string& uri, OrcaData& );

private:
    std::string arrangement_{ "T" };
};

}  // namespace atlas::orca


//------------------------------------------------------------------------------------------------------
