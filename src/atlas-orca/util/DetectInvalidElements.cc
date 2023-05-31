/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DetectInvalidElements.h"

#include <cmath>

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {
namespace orca {

bool DetectInvalidElement::invalid_quad_2d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                        const PointLonLat& p_NW ) const {
    double dlat_W = p_NW.lat() - p_SW.lat();
    double dlat_E = p_NE.lat() - p_SE.lat();
    double dlon_N = p_NE.lon() - p_NW.lon();
    double dlon_S = p_SE.lon() - p_SW.lon();
    return ( dlat_W < -1.e-10 || dlat_E < -1.e-10 || dlon_N < -1.e-10 || dlon_S < -1.e-10 );
}

bool DetectInvalidElement::invalid_quad_3d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                        const PointLonLat& p_NW ) const {
    PointXYZ xyz_SW{sphere_.xyz( p_SW )};
    PointXYZ xyz_SE{sphere_.xyz( p_SE )};
    PointXYZ xyz_NE{sphere_.xyz( p_NE )};
    PointXYZ xyz_NW{sphere_.xyz( p_NW )};
    return not atlas::interpolation::element::Quad3D{xyz_SW, xyz_SE, xyz_NE, xyz_NW}.validate();
}

bool DetectInvalidElement::diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                            const PointLonLat& p_NW, double largest_diagonal ) const {
    if ( largest_diagonal == 0. ) {
        return false;
    }
    // example with extended grids where grid folds over itself, or ORCA2 grids
    double diagonal_treshold_2 = largest_diagonal * largest_diagonal;
    double d2_NW_SE            = PointLonLat::distance2( p_NW, p_SE );
    double d2_SW_NE            = PointLonLat::distance2( p_SW, p_NE );
    return std::max( d2_NW_SE, d2_SW_NE ) > diagonal_treshold_2;
}

bool DetectInvalidElement::diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                            const PointLonLat& p_NW ) const {
    return diagonal_too_large( p_SW, p_SE, p_NE, p_NW, largest_diatonal_ );
}

bool DetectInvalidElement::invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                        const PointLonLat& p_NW, Statistics& statistics ) const {
    double lat_max = std::max( {p_SW.lat(), p_SE.lat(), p_NE.lat(), p_NW.lat()} );

    if ( invalid_quad_2d( p_SW, p_SE, p_NE, p_NW ) && lat_max < 45. ) {
        statistics.invalid_quads_2d++;
        statistics.invalid_elements++;
        return true;
    }
    if ( lat_max < 60. && diagonal_too_large( p_SW, p_SE, p_NE, p_NW ) ) {
        statistics.diagonal_too_large++;
        statistics.invalid_elements++;
        return true;
    }
    if ( invalid_quad_3d( p_SW, p_SE, p_NE, p_NW ) ) {
        statistics.invalid_quads_3d++;
        statistics.invalid_elements++;
        return true;
    }
    if ( orca2_ ) {
        if ( lat_max < 60. ) {
            constexpr util::NormaliseLongitude normalized{-180.};
            double lon_min = normalized( std::min( {p_SW.lon(), p_SE.lon(), p_NE.lon(), p_NW.lon()} ) );
            if ( lon_min > -20. && lon_min < 20. ) {
                double dlat_W = p_NW.lat() - p_SW.lat();
                double dlat_E = p_NE.lat() - p_SE.lat();
                double dlon_N = p_NE.lon() - p_NW.lon();
                double dlon_S = p_SE.lon() - p_SW.lon();
                if ( std::max( dlon_N, dlon_S ) > 2. * std::min( dlon_N, dlon_S ) ) {
                    statistics.invalid_elements++;
                    return true;
                }
            }
        }
    }
    return false;
}

bool DetectInvalidElement::invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                        const PointLonLat& p_NW ) const {
    Statistics statistics;
    return invalid_element( p_SW, p_SE, p_NE, p_NW, statistics );
}


}  // namespace orca
}  // namespace atlas
