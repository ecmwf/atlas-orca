/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "OrcaData.h"

#if ATLAS_ORCA_HAVE_ECKIT_CODEC
#include "eckit/codec/codec.h"
#else
// Backward compatibility, DEPRECATED!
#include "atlas/io/atlas-io.h"
namespace eckit::codec {
using RecordWriter = atlas::io::RecordWriter;
using atlas::io::ArrayReference;
using atlas::io::ref;
}  // namespace eckit::codec
#endif

#include "eckit/log/ProgressTimer.h"

#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "atlas-orca/util/ComputeUid.h"
#include "atlas-orca/util/DetectInvalidElements.h"
#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"
#include "atlas-orca/util/OrcaPeriodicity.h"


namespace atlas::orca {

void OrcaData::setGhost() {
    idx_t ni = dimensions[0];
    idx_t nj = dimensions[1];

    OrcaPeriodicity compute_master{ *this };
    for ( idx_t j = 0; j < nj; ++j ) {
        for ( idx_t i = 0; i < ni; ++i ) {
            idx_t n = ni * j + i;
            Flag flag{ flags[n] };
            flag.unset( Flag::GHOST );
            auto master = compute_master( i, j );
            if ( j < halo[HALO_SOUTH] ) {
                flag.set( Flag::GHOST );
            }
            else if ( master.i != i || master.j != j ) {
                flag.set( Flag::GHOST );
            }
        }
    }
}

void OrcaData::makeHaloConsistent() {
    size_t ni = dimensions[0];
    size_t nj = dimensions[1];
    OrcaPeriodicity compute_master{ *this };
    for ( size_t j = 0; j < nj; ++j ) {
        for ( size_t i = 0; i < ni; ++i ) {
            size_t n        = ni * j + i;
            auto master     = compute_master( static_cast<idx_t>( i ), static_cast<idx_t>( j ) );
            size_t n_master = ni * master.j + master.i;
            if ( n_master != n ) {
                if ( lon[n] != lon[n_master] || lat[n] != lat[n_master] ) {
                    PointLonLat point_n{ lon[n], lat[n] };
                    PointLonLat point_n_master{ lon[n_master], lat[n_master] };
                    double distance = geometry::Earth().distance( point_n, point_n_master );
                    if ( distance > 1 /*metre*/ ) {
                        Log::warning() << "Fixed halo inconsistency for {" << i << "," << j << "}:  " << point_n
                                       << " --> " << point_n_master << "   ( distance = " << distance << " m )"
                                       << std::endl;
                    }
                }
                lon[n] = lon[n_master];
                lat[n] = lat[n_master];

                Flag flag_n{ flags[n] };
                Flag flag_n_master{ flags[n_master] };
                if ( flag_n.test( Flag::WATER ) != flag_n_master.test( Flag::WATER ) ) {
                    flag_n.unset( Flag::WATER );
                    if ( flag_n_master.test( Flag::WATER ) ) {
                        flag_n.set( Flag::WATER );
                    }
                    Log::warning() << "Fixed halo inconsistency for {" << i << "," << j << "}: land mask" << std::endl;
                }
            }
        }
    }
}

std::string atlas::orca::OrcaData::computeUid( const util::Config& config ) {
    checkSetup();

    std::string arrangement;
    if ( not config.get( "arrangement", arrangement ) ) {
        config.get( "orca_arrangement", arrangement );
    }
    size_t size = lon.size();
    return orca::compute_uid( arrangement, lon.data(), lat.data(), size );
}

DetectInvalidElement::Statistics atlas::orca::OrcaData::detectInvalidElements( const util::Config& config ) {
    checkSetup();

    DetectInvalidElement::Statistics stats;

    ATLAS_TRACE( "Detect invalid elements" );
    double diagonal_factor = config.getDouble( "diagonal-factor", 2.5 );

    idx_t ni          = dimensions[0];
    idx_t nj          = dimensions[1];
    idx_t nx          = ni - halo[HALO_EAST] - halo[HALO_WEST];
    double resolution = 360. / static_cast<double>( nx );

    util::Config detection_config;
    detection_config.set( "ORCA2", ( std::abs( resolution - 2. ) < 0.1 ) );
    detection_config.set( "diagonal", resolution * diagonal_factor );
    DetectInvalidElement detect( detection_config );

    auto is_water = [&]( int n ) -> bool { return Flag{ flags[n] }.test( Flag::WATER ); };

    eckit::ProgressTimer progress( "Detect invalid elements", ( nj - 1 ) * ( ni - 1 ), "point", 5., Log::trace() );
    for ( idx_t j = 0; j < nj - 1; ++j ) {
        for ( idx_t i = 0; i < ni - 1; ++i ) {
            ++progress;
            idx_t n_SW = j * ni + i;
            idx_t n_SE = n_SW + 1;
            idx_t n_NW = ( j + 1 ) * ni + i;
            idx_t n_NE = n_NW + 1;

            Flag{ flags[n_SW] }.unset( Flag::INVALID_ELEMENT );

            bool element_contains_water = is_water( n_SW ) || is_water( n_SE ) || is_water( n_NW ) || is_water( n_NE );

            PointLonLat p_SW{ lon[n_SW], lat[n_SW] };
            PointLonLat p_SE{ lon[n_SE], lat[n_SE] };
            PointLonLat p_NW{ lon[n_NW], lat[n_NW] };
            PointLonLat p_NE{ lon[n_NE], lat[n_NE] };
            p_SE.normalise( p_SW.lon() - 180. );
            p_NW.normalise( p_SW.lon() - 180. );
            p_NE.normalise( p_SW.lon() - 180. );

            if ( detect.invalid_element( p_SW, p_SE, p_NE, p_NW, stats ) ) {
                if ( element_contains_water ) {
                    // So far this only occurs for ORCA2_F grid
                    int invalid_factor = 10;
                    if ( detect.diagonal_too_large( p_SW, p_SE, p_NE, p_NW, invalid_factor * resolution ) ) {
                        Log::warning() << "Element {I,J} = {" << i << "," << j
                                       << "} is invalidated although it contains water!" << std::endl;
                        Log::warning() << "  South-West point: {lon,lat} = " << p_SW << std::endl;
                        Log::warning() << "  diagonal > " << invalid_factor << " * ref_length" << std::endl;
                    }
                    else {
                        for ( int k = invalid_factor - 1; k >= diagonal_factor; --k ) {
                            if ( detect.diagonal_too_large( p_SW, p_SE, p_NE, p_NW,
                                                            static_cast<double>( k ) * resolution ) ) {
                                Log::warning() << "Element {I,J} = {" << i << "," << j
                                               << "} is not invalidated as it contains water, even though (diagonal > "
                                               << diagonal_factor << " * ref_length)." << std::endl;
                                Log::warning() << "  South-West point: {lon,lat} = " << p_SW << std::endl;
                                Log::warning() << "  diagonal > " << k << " * ref_length" << std::endl;
                                break;
                            }
                        }
                        continue;
                    }
                }
                Flag{ flags[n_SW] }.set( Flag::INVALID_ELEMENT );
            }
        }
    }
    return stats;
}

void atlas::orca::OrcaData::checkSetup() {
    ATLAS_ASSERT( dimensions[0] >= 0 );
    ATLAS_ASSERT( dimensions[1] >= 0 );
    ATLAS_ASSERT( halo[0] >= 0 );
    ATLAS_ASSERT( halo[1] >= 0 );
    ATLAS_ASSERT( halo[2] >= 0 );
    ATLAS_ASSERT( halo[3] >= 0 );
    ATLAS_ASSERT( pivot[0] >= 0 );
    ATLAS_ASSERT( pivot[1] >= 0 );
    ATLAS_ASSERT( not lon.empty() );
    ATLAS_ASSERT( not lat.empty() );
    ATLAS_ASSERT( not flags.empty() );
    size_t size = dimensions[0] * dimensions[1];
    ATLAS_ASSERT( lon.size() == size );
    ATLAS_ASSERT( lat.size() == size );
    ATLAS_ASSERT( flags.size() == size );
}

size_t atlas::orca::OrcaData::write( const eckit::PathName& path, const util::Config& config ) {
    checkSetup();
    eckit::codec::RecordWriter record;
    record.compression( config.getString( "compression", "none" ) );
    record.set( "version", 0 );
    record.set( "dimensions", eckit::codec::ref( dimensions ) );
    record.set( "halo", eckit::codec::ref( halo ) );
    record.set( "pivot", eckit::codec::ref( pivot ) );
    record.set( "longitude", eckit::codec::ArrayReference( lon.data(), dimensions ) );
    record.set( "latitude", eckit::codec::ArrayReference( lat.data(), dimensions ) );
    record.set( "flags", eckit::codec::ArrayReference( flags.data(), dimensions ) );
    return record.write( path );
}


}  // namespace atlas::orca
