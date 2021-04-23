/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "OrcaData.h"

#include "atlas/io/atlas-io.h"

#include "atlas-orca/util/OrcaPeriodicity.h"
#include "atlas-orca/util/DetectInvalidElements.h"
#include "atlas-orca/util/ComputeUid.h"
#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"

namespace atlas {
namespace orca {

void OrcaData::setGhost() {
    idx_t ni = dimensions[0];
    idx_t nj = dimensions[1];

    OrcaPeriodicity compute_master{*this};
    for( idx_t j = 0; j < nj; ++j ) {
        for( idx_t i = 0; i < ni; ++i ) {
            idx_t n = ni*j + i;
            Flag flag{flags[n]};
            flag.unset(Flag::GHOST);
            auto master = compute_master(i,j);
            if( j < halo[HALO_SOUTH] ) {
                flag.set(Flag::GHOST);
            }
            else if( master.i != i || master.j != j ) {
                flag.set(Flag::GHOST);
            }
        }
    }
}

std::string atlas::orca::OrcaData::computeUid(const util::Config &config) {
    checkSetup();

    std::string arrangement;
    if( not config.get("arrangement",arrangement) ) {
        config.get("orca_arrangement",arrangement);
    }
    return orca::compute_uid(arrangement,lon,lat,dimensions[0],dimensions[1],halo);
}

size_t atlas::orca::OrcaData::detectInvalidElementss(const util::Config &config) {
    checkSetup();

    size_t invalid_elements = 0;

    ATLAS_TRACE("Detect invalid elements");
    double diagonal_factor = config.getDouble("diagonal-factor",2.5);

    idx_t ni = dimensions[0];
    idx_t nj = dimensions[1];
    idx_t nx = ni - halo[HALO_EAST] - halo[HALO_WEST];
    double resolution = 360. / double(nx);

    util::Config detection_config;
    detection_config.set("ORCA2", (std::abs(resolution-2.) < 0.1) );
    detection_config.set("diagonal", resolution*diagonal_factor );
    atlas::orca::meshgenerator::DetectInvalidElement detect(detection_config);

    auto is_water = [&](int n) -> bool {
        return Flag{flags[n]}.test(Flag::WATER);
    };

    for( idx_t j=0; j<nj-1; ++j ) {
        for( idx_t i=0; i<ni-1; ++i ) {
            idx_t n_SW = j*ni + i;
            idx_t n_SE = n_SW + 1;
            idx_t n_NW = (j+1)*ni + i;
            idx_t n_NE = n_NW + 1;

            Flag{flags[n_SW]}.unset(Flag::INVALID_ELEMENT);

            bool element_contains_water = is_water(n_SW) || is_water(n_SE) || is_water(n_NW) || is_water(n_NE);

            PointLonLat p_SW{ lon[n_SW], lat[n_SW] };
            PointLonLat p_SE{ lon[n_SE], lat[n_SE] };
            PointLonLat p_NW{ lon[n_NW], lat[n_NW] };
            PointLonLat p_NE{ lon[n_NE], lat[n_NE] };
            p_SE.normalise(p_SW.lon()-180.);
            p_NW.normalise(p_SW.lon()-180.);
            p_NE.normalise(p_SW.lon()-180.);

            if(  detect.invalid_element(p_SW,p_SE,p_NE,p_NW) ) {
                if( element_contains_water ) {
                    // So far this only occurs for ORCA2_F grid
                    int invalid_factor = 10;
                    if( detect.diagonal_too_large(p_SW,p_SE,p_NE,p_NW,invalid_factor*resolution) ) {
                        Log::warning() << "Element {I,J} = {"<<i<<","<<j<<"} is invalidated although it contains water!" << std::endl;
                        Log::warning() << "  South-West point: {lon,lat} = " << p_SW << std::endl;
                        Log::warning() << "  diagonal > " << invalid_factor << " * ref_length" << std::endl;
                    }
                    else {
                        for( int k=invalid_factor-1; k>=diagonal_factor; --k ) {
                            if( detect.diagonal_too_large(p_SW,p_SE,p_NE,p_NW,double(k)*resolution) ) {
                                Log::warning() << "Element {I,J} = {"<<i<<","<<j<<"} is not invalidated as it contains water, even though (diagonal > " << diagonal_factor << " * ref_length)." << std::endl;
                                Log::warning() << "  South-West point: {lon,lat} = " << p_SW << std::endl;
                                Log::warning() << "  diagonal > " << k << " * ref_length" << std::endl;
                                break;
                            }
                        }
                        continue;
                    }
                }
                Flag{flags[n_SW]}.set(Flag::INVALID_ELEMENT);
                ++invalid_elements;
            }
        }
    }
    return invalid_elements;
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
}

size_t atlas::orca::OrcaData::write(const eckit::PathName &path, const util::Config &config) {
    checkSetup();
    ATLAS_ASSERT( lon.size() == dimensions[0]*dimensions[1] );
    ATLAS_ASSERT( lat.size() == dimensions[0]*dimensions[1] );
    ATLAS_ASSERT( flags.size() == dimensions[0]*dimensions[1]);
    io::RecordWriter record;
    record.compression( config.getString("compression","none") );
    record.set("version",0);
    record.set("dimensions",io::ref(dimensions));
    record.set("halo", io::ref(halo));
    record.set("pivot", io::ref(pivot) );
    record.set("longitude", io::ArrayReference(lon.data(), dimensions));
    record.set("latitude", io::ArrayReference(lat.data(), dimensions));
    record.set("flags", io::ArrayReference(flags.data(), dimensions));
    return record.write(path);
}


}  // namespace orca
}  // namespace atlas
