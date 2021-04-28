/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ComputeUid.h"

#include <algorithm>

#include "eckit/utils/ByteSwap.h"
#include "eckit/utils/MD5.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/vector.h"

#include "atlas-orca/util/Enums.h"

namespace atlas {
namespace orca {

namespace uid_v0 {
std::string compute_uid( const std::string& arrangement, const double lon[], const double lat[], int nx_halo,
                         int ny_halo, const int32_t halo[] ) {
    ATLAS_TRACE();

    bool valid_arrangement = false;
    if ( arrangement == "F" ) {
        valid_arrangement = true;
    }
    if ( arrangement == "T" ) {
        valid_arrangement = true;
    }
    if ( arrangement == "U" ) {
        valid_arrangement = true;
    }
    if ( arrangement == "V" ) {
        valid_arrangement = true;
    }
    if ( arrangement == "W" ) {
        valid_arrangement = true;
    }
    if ( arrangement.empty() ) {
        valid_arrangement = true;
    }
    ATLAS_ASSERT( valid_arrangement );

    idx_t i_begin = halo[HALO_WEST];
    idx_t i_end   = nx_halo - halo[HALO_EAST];
    idx_t j_begin = halo[HALO_SOUTH];
    idx_t j_end   = ny_halo - halo[HALO_NORTH];
    idx_t nx      = i_end - i_begin;
    idx_t ny      = j_end - j_begin;
    atlas::vector<std::int32_t> lonlat( 2 * nx * ny );
    size_t c{0};
    for ( idx_t j = j_begin; j < j_end; ++j ) {
        for ( idx_t i = i_begin; i < i_end; ++i ) {
            idx_t n     = j * nx_halo + i;
            lonlat[c++] = util::microdeg( lon[n] );
            lonlat[c++] = util::microdeg( lat[n] );
        }
    }
    eckit::MD5 hasher;

    hasher.add( lonlat.data(), lonlat.size() * long( sizeof( std::int32_t ) ) );

    hasher.add( halo, 4 * sizeof( std::int32_t ) );

    if ( arrangement == "W" ) {
        hasher.add( arrangement );
    }
    return hasher.digest();
}

std::string compute_uid( const std::string& arrangement, const std::vector<double>& lon, const std::vector<double>& lat,
                         int nx_halo, int ny_halo, const std::array<std::int32_t, 4>& halo ) {
    return compute_uid( arrangement, lon.data(), lat.data(), nx_halo, ny_halo, halo.data() );
}
}  // namespace uid_v0

std::string compute_uid( const std::string& arrangement, const double lon[], const double lat[], int nx_halo,
                         int ny_halo, const int32_t halo[] ) {
    ATLAS_TRACE();

    std::string P = arrangement;
    if ( P != "T" && P != "W" && P != "F" && P != "U" && P != "V" ) {
        ATLAS_THROW_EXCEPTION( "arrangement expected to be any of {F,T,U,V,W}. Received: " << P );
    }

    eckit::MD5 hasher;
    hasher.add( P );

    idx_t size = nx_halo * ny_halo;

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

std::string compute_uid( const std::string& arrangement, const std::vector<double>& lon, const std::vector<double>& lat,
                         int nx_halo, int ny_halo, const std::array<std::int32_t, 4>& halo ) {
    return compute_uid( arrangement, lon.data(), lat.data(), nx_halo, ny_halo, halo.data() );
}


}  // namespace orca
}  // namespace atlas
