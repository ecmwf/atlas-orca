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

#include "atlas/util/Config.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace orca {

class DetectInvalidElement {
public:
    struct Statistics {
        size_t invalid_elements{0};
        size_t invalid_quads_3d{0};
        size_t invalid_quads_2d{0};
        size_t diagonal_too_large{0};
    };

    DetectInvalidElement( const util::Config& config ) {
        config.get( "ORCA2", orca2_ );
        config.get( "diagonal", largest_diatonal_ );
    }

    bool invalid_quad_2d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

    bool invalid_quad_3d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

    bool diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                             const PointLonLat& p_NW, double largest_diagonal ) const;

    bool diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                             const PointLonLat& p_NW ) const;

    bool invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW, Statistics& statistics ) const;

    bool invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

private:
    geometry::Earth sphere_;
    double largest_diatonal_{0};
    bool orca2_{false};
};

}  // namespace orca
}  // namespace atlas
