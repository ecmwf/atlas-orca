/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ComputeUid.h"

#include <algorithm>

#include "eckit/eckit_config.h"
#include "eckit/utils/ByteSwap.h"
#include "eckit/utils/MD5.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/vector.h"

#include "atlas-orca/util/Enums.h"

namespace atlas {
namespace orca {

std::string compute_uid( const std::string& arrangement, const double lon[], const double lat[], size_t size ) {
    ATLAS_TRACE();

    std::string P = arrangement;
    if ( P != "T" && P != "W" && P != "F" && P != "U" && P != "V" ) {
        ATLAS_THROW_EXCEPTION( "arrangement expected to be any of {F,T,U,V,W}. Received: " << P );
    }

    eckit::MD5 hasher;
    hasher.add( P );

    if ( eckit_LITTLE_ENDIAN ) {
        hasher.add( lat, size * sizeof( double ) );
        hasher.add( lon, size * sizeof( double ) );
    }
    else {
        atlas::vector<double> latitude( size );
        atlas::vector<double> longitude( size );
        for ( idx_t n = 0; n < size; ++n ) {
            latitude[n]  = lat[n];
            longitude[n] = lon[n];
        }
        eckit::byteswap( latitude.data(), size );
        eckit::byteswap( longitude.data(), size );
        hasher.add( latitude.data(), size * sizeof( double ) );
        hasher.add( longitude.data(), size * sizeof( double ) );
    }
    return hasher.digest();
}

}  // namespace orca
}  // namespace atlas
