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


namespace atlas::orca {

class DetectInvalidElement {
public:
    enum class Reason {
        none,
        invalid_quad_2d,
        invalid_quad_3d,
        diagonal_too_large,
        zero_diagonal,
        zero_edge,
        self_intersecting,
        invalid_jacobian,
        poor_quality,
        southern_edge_aspect,
        eorca025_southern_cap,
        orca2_longitude_aspect,
        western_europe_edge_aspect
    };

    struct Statistics {
        Reason last_reason{ Reason::none };
        size_t invalid_elements{ 0 };
        size_t invalid_quads_3d{ 0 };
        size_t invalid_quads_2d{ 0 };
        size_t diagonal_too_large{ 0 };
        size_t zero_diagonal{ 0 };
        size_t zero_edge{ 0 };
        size_t self_intersecting{ 0 };
        size_t invalid_jacobian{ 0 };
        size_t poor_quality{ 0 };
        size_t southern_edge_aspect{ 0 };
        size_t western_europe_edge_aspect{ 0 };
        double average_diagonal{ 0 };
        double minimum_diagonal{ 1.e30 };
        double maximum_diagonal{ 0 };
        size_t num_diagonals{ 0 };
    };

    static const char* reasonString( Reason reason );

    explicit DetectInvalidElement( const util::Config& config );

    bool invalid_quad_2d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

    bool invalid_quad_3d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

    bool invalid_quad_3d( const PointXYZ& p_SW, const PointXYZ& p_SE, const PointXYZ& p_NE,
                          const PointXYZ& p_NW ) const;

    bool diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                             const PointLonLat& p_NW, double largest_diagonal ) const;

    bool diagonal_too_large( const PointXYZ& p_SW, const PointXYZ& p_SE, const PointXYZ& p_NE,
                             const PointXYZ& p_NW, double largest_diagonal ) const;

    bool diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                             const PointLonLat& p_NW ) const;

    bool invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW, Statistics& statistics ) const;

    bool invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                          const PointLonLat& p_NW ) const;

private:
    geometry::Earth sphere_;
    double length_tolerance_;
    double largest_diagonal_{ 0 };
    bool ORCA2_{ false };
    bool eORCA025_{ false };
};

}  // namespace atlas::orca
