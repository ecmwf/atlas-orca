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

#include <string>

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/URI.h"

namespace atlas::orca {

class OrcaDataFile {
public:
    explicit OrcaDataFile( const std::string& uri, const std::string& checksum = "" );

    const eckit::PathName& path() const { return path_; }

    explicit operator const eckit::PathName&() const { return path(); }
    operator std::string() const { return path_.asString(); }
    explicit operator const char*() const { return c_str(); }
    const char* c_str() const { return path_.localPath(); }

private:
    eckit::URI uri_;
    eckit::PathName path_;
    std::string checksum_;
};

}  // namespace atlas::orca
