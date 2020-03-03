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


#include "atlas-data/Library.h"


namespace atlas {
namespace data {


REGISTER_LIBRARY( Library );


Library::Library() : eckit::system::Library( "atlas-data" ) {}


const Library& Library::instance() {
    static Library library;
    return library;
}


const void* Library::addr() const {
    return this;
}


std::string Library::version() const {
    return "not available";
}


std::string Library::gitsha1( unsigned int ) const {
    return "not available";
}


}  // namespace data
}  // namespace atlas
