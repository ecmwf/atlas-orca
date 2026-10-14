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

#include <algorithm>
#include <array>
#include <cmath>

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/util/NormaliseLongitude.h"


namespace atlas::orca {

const char* DetectInvalidElement::reasonString( Reason reason ) {
    switch ( reason ) {
        case Reason::none: return "none";
        case Reason::invalid_quad_2d: return "invalid quad in longitude-latitude space";
        case Reason::invalid_quad_3d: return "invalid quad in 3D space";
        case Reason::diagonal_too_large: return "diagonal too large";
        case Reason::zero_diagonal: return "zero-length diagonal";
        case Reason::zero_edge: return "zero-length edge";
        case Reason::self_intersecting: return "self-intersecting";
        case Reason::invalid_jacobian: return "invalid Jacobian";
        case Reason::poor_quality: return "poor quality";
        case Reason::southern_edge_aspect: return "southern edge aspect ratio";
        case Reason::eorca025_southern_cap: return "eORCA025 southern cap";
        case Reason::orca2_longitude_aspect: return "ORCA2 longitude aspect ratio";
        case Reason::western_europe_edge_aspect: return "western Europe edge aspect ratio";
    }
    return "unknown";
}

namespace {

/**
 * @brief Converts an angular degree segment on a sphere to its straight-line 3D Euclidean distance.
 * @param angular_separation_degrees The total angle between two points (0.0 to 180.0).
 * @param R The radius of the sphere in meters.
 */
inline double degree_segment_to_euclidean(double angular_separation_degrees, double R) noexcept {
    // PI / 360 handles both dividing the angle by 2 and converting degrees to radians
    constexpr double deg_to_half_rad = M_PI / 360.0;

    return 2.0 * R * std::sin(angular_separation_degrees * deg_to_half_rad);
}

/**
 * @brief Converts a 3D straight-line Euclidean distance back into an angular degree segment on a sphere.
 * @param euclidean_distance The physical straight-line distance between the points.
 * @param R The radius of the sphere.
 */
inline double euclidean_to_degree_segment(double euclidean_distance, double R) noexcept {
    // Clamping the input ratio between -1.0 and 1.0 prevents NaN floating-point
    // crashes if slight numerical precision leakage makes (dist / 2R) slightly exceed 1.0.
    double sin_half_angle = std::clamp(euclidean_distance / (2.0 * R), -1.0, 1.0);

    // 360 / PI handles both doubling the half-angle and converting radians to degrees
    constexpr double rad_to_double_deg = 360.0 / M_PI;

    return std::asin(sin_half_angle) * rad_to_double_deg;
}

bool is_zero_length(const Point3& p1, const Point3& p2, double epsilon = 1.e-14) noexcept {
    // 1. Compute individual component deltas
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    // 2. The dot product (A-C) · (A-C)
    const double squared_distance = (dx * dx) + (dy * dy) + (dz * dz);

    // 3. Evaluate against a strict floating-point tolerance near zero
    // For standard double-precision, 1e-14 is extremely safe and precise.
    // constexpr double epsilon = 1e-14;

    return squared_distance < epsilon*epsilon;
};

double diagonal_in_degrees(const Point3& p1, const Point3& p2, double radius) noexcept {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];
    double diagonal_2 = (dx * dx) + (dy * dy) + (dz * dz);
    double diagonal = std::sqrt(diagonal_2);
    return euclidean_to_degree_segment(diagonal, radius);
}

template <typename Point>
bool has_zero_diagonal(const Point& p1, const Point& p2, const Point& p3, const Point& p4, double epsilon = 1.e-14) noexcept {
    return (is_zero_length(p1, p3, epsilon) || is_zero_length(p2, p4, epsilon));
}

template <typename Point>
[[maybe_unused]] bool has_zero_edge( const Point& p1, const Point& p2, const Point& p3, const Point& p4,
                    double epsilon = 1.e-14 ) noexcept {
    return is_zero_length( p1, p2, epsilon ) || is_zero_length( p2, p3, epsilon ) ||
           is_zero_length( p3, p4, epsilon ) || is_zero_length( p4, p1, epsilon );
}

bool has_self_intersection( const Point3& p_SW, const Point3& p_SE, const Point3& p_NE,
                            const Point3& p_NW ) noexcept {
    Point3 radial{ p_SW[0] + p_SE[0] + p_NE[0] + p_NW[0], p_SW[1] + p_SE[1] + p_NE[1] + p_NW[1],
                   p_SW[2] + p_SE[2] + p_NE[2] + p_NW[2] };
    const double radial_length =
        std::sqrt( radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2] );
    if ( radial_length == 0. ) {
        return true;
    }
    for ( size_t component = 0; component < 3; ++component ) {
        radial[component] /= radial_length;
    }

    const Point3 reference = std::abs( radial[2] ) < 0.9 ? Point3{ 0., 0., 1. } : Point3{ 1., 0., 0. };
    Point3 axis_x{ reference[1] * radial[2] - reference[2] * radial[1],
                   reference[2] * radial[0] - reference[0] * radial[2],
                   reference[0] * radial[1] - reference[1] * radial[0] };
    const double axis_x_length =
        std::sqrt( axis_x[0] * axis_x[0] + axis_x[1] * axis_x[1] + axis_x[2] * axis_x[2] );
    for ( size_t component = 0; component < 3; ++component ) {
        axis_x[component] /= axis_x_length;
    }
    const Point3 axis_y{ radial[1] * axis_x[2] - radial[2] * axis_x[1],
                         radial[2] * axis_x[0] - radial[0] * axis_x[2],
                         radial[0] * axis_x[1] - radial[1] * axis_x[0] };

    auto project = [&]( const Point3& point ) {
        return Point2{ point[0] * axis_x[0] + point[1] * axis_x[1] + point[2] * axis_x[2],
                       point[0] * axis_y[0] + point[1] * axis_y[1] + point[2] * axis_y[2] };
    };
    auto orientation = []( const Point2& a, const Point2& b, const Point2& c ) {
        return ( b[0] - a[0] ) * ( c[1] - a[1] ) - ( b[1] - a[1] ) * ( c[0] - a[0] );
    };
    auto segments_intersect = [&]( const Point2& a, const Point2& b, const Point2& c, const Point2& d ) {
        constexpr double tolerance = 1.e-12;
        const double abc = orientation( a, b, c );
        const double abd = orientation( a, b, d );
        const double cda = orientation( c, d, a );
        const double cdb = orientation( c, d, b );
        return abc * abd < -tolerance && cda * cdb < -tolerance;
    };

    const Point2 sw = project( p_SW );
    const Point2 se = project( p_SE );
    const Point2 ne = project( p_NE );
    const Point2 nw = project( p_NW );
    return segments_intersect( sw, se, ne, nw ) || segments_intersect( se, ne, nw, sw );
}

double scaled_jacobian( const Point3& origin, const Point3& edge_i, const Point3& edge_j ) noexcept {
    // Measures the signed corner quality relative to the outward radial direction:
    //     ((edge_i - origin) x (edge_j - origin)) . origin
    //     --------------------------------------------------
    //       |edge_i - origin| |edge_j - origin| |origin|
    //
    // Its magnitude is the sine of the angle between the adjacent edges after their
    // cross product is projected onto the local radial direction. A value of 1 is a
    // right-angled corner with the expected orientation, 0 means that the edges are
    // collinear or that their corner normal is tangent to the sphere, and a negative
    // value indicates reversed local orientation.

    const Point3 tangent_i{ edge_i[0] - origin[0], edge_i[1] - origin[1], edge_i[2] - origin[2] };
    const Point3 tangent_j{ edge_j[0] - origin[0], edge_j[1] - origin[1], edge_j[2] - origin[2] };
    const Point3 normal{ tangent_i[1] * tangent_j[2] - tangent_i[2] * tangent_j[1],
                         tangent_i[2] * tangent_j[0] - tangent_i[0] * tangent_j[2],
                         tangent_i[0] * tangent_j[1] - tangent_i[1] * tangent_j[0] };
    const double tangent_i_length = std::sqrt( tangent_i[0] * tangent_i[0] + tangent_i[1] * tangent_i[1] +
                                               tangent_i[2] * tangent_i[2] );
    const double tangent_j_length = std::sqrt( tangent_j[0] * tangent_j[0] + tangent_j[1] * tangent_j[1] +
                                               tangent_j[2] * tangent_j[2] );
    const double origin_length =
        std::sqrt( origin[0] * origin[0] + origin[1] * origin[1] + origin[2] * origin[2] );
    const double scale = tangent_i_length * tangent_j_length * origin_length;
    return scale == 0. ? 0. : ( normal[0] * origin[0] + normal[1] * origin[1] + normal[2] * origin[2] ) / scale;
}

std::array<double, 4> scaled_jacobians( const Point3& p_SW, const Point3& p_SE, const Point3& p_NE,
                                        const Point3& p_NW ) noexcept {
    return { scaled_jacobian( p_SW, p_SE, p_NW ), scaled_jacobian( p_SE, p_NE, p_SW ),
             scaled_jacobian( p_NE, p_NW, p_SE ), scaled_jacobian( p_NW, p_SW, p_NE ) };
}

[[maybe_unused]]bool has_invalid_jacobian( const Point3& p_SW, const Point3& p_SE, const Point3& p_NE,
                           const Point3& p_NW ) noexcept {
    constexpr double tolerance = 1.e-8;
    const auto jacobians = scaled_jacobians( p_SW, p_SE, p_NE, p_NW );

    // A valid non-degenerate quad has a consistent local orientation at every corner.
    // A Jacobian close to zero identifies a collapsed or locally tangent corner. Mixed
    // signs mean that the orientation changes within the quad, indicating a concavity,
    // fold, or self-intersection. Four negative values are consistent reversed winding
    // and are deliberately not classified as invalid by this check alone.
    bool positive = false;
    bool negative = false;
    for ( double jacobian : jacobians ) {
        if ( std::abs( jacobian ) <= tolerance ) {
            return true;
        }
        positive = positive || jacobian > 0.;
        negative = negative || jacobian < 0.;
    }
    return positive && negative;
}

double chord_length( const Point3& a, const Point3& b ) noexcept {
    const double dx = b[0] - a[0];
    const double dy = b[1] - a[1];
    const double dz = b[2] - a[2];
    return std::sqrt( dx * dx + dy * dy + dz * dz );
}

double element_quality( const Point3& p_SW, const Point3& p_SE, const Point3& p_NE,
                        const Point3& p_NW ) noexcept {
    const auto jacobians = scaled_jacobians( p_SW, p_SE, p_NE, p_NW );
    double corner_quality = 1.;
    for ( double jacobian : jacobians ) {
        corner_quality = std::min( corner_quality, std::abs( jacobian ) );
    }

    const std::array<double, 4> edge_lengths{ chord_length( p_SW, p_SE ), chord_length( p_SE, p_NE ),
                                               chord_length( p_NE, p_NW ), chord_length( p_NW, p_SW ) };
    const auto edge_extrema = std::minmax_element( edge_lengths.begin(), edge_lengths.end() );
    const double edge_quality = *edge_extrema.second == 0. ? 0. : *edge_extrema.first / *edge_extrema.second;
    return std::min( corner_quality, edge_quality );
}

}

DetectInvalidElement::DetectInvalidElement( const util::Config& config ) :
    length_tolerance_( degree_segment_to_euclidean( 1.e-4, sphere_.radius() ) ) {
    config.get( "ORCA2", ORCA2_ );
    config.get( "eORCA025", eORCA025_ );
    config.get( "diagonal", largest_diagonal_ );
}


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
    const PointXYZ xyz_SW{ sphere_.xyz( p_SW ) };
    const PointXYZ xyz_SE{ sphere_.xyz( p_SE ) };
    const PointXYZ xyz_NE{ sphere_.xyz( p_NE ) };
    const PointXYZ xyz_NW{ sphere_.xyz( p_NW ) };
    return invalid_quad_3d( xyz_SW, xyz_SE, xyz_NE, xyz_NW );
}

bool DetectInvalidElement::invalid_quad_3d( const PointXYZ& p_SW, const PointXYZ& p_SE, const PointXYZ& p_NE,
                                            const PointXYZ& p_NW ) const {
    return not atlas::interpolation::element::Quad3D{ p_SW, p_SE, p_NE, p_NW }.validate();
}

bool DetectInvalidElement::diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE,
                                               const PointLonLat& p_NE, const PointLonLat& p_NW,
                                               double largest_diagonal ) const {
    const PointXYZ xyz_SW{ sphere_.xyz( p_SW ) };
    const PointXYZ xyz_SE{ sphere_.xyz( p_SE ) };
    const PointXYZ xyz_NE{ sphere_.xyz( p_NE ) };
    const PointXYZ xyz_NW{ sphere_.xyz( p_NW ) };
    return diagonal_too_large( xyz_SW, xyz_SE, xyz_NE, xyz_NW, largest_diagonal );
}

bool DetectInvalidElement::diagonal_too_large( const PointXYZ& p_SW, const PointXYZ& p_SE,
                                               const PointXYZ& p_NE, const PointXYZ& p_NW,
                                               double largest_diagonal ) const {
    if ( largest_diagonal == 0. ) {
        return false;
    }
    // example with extended grids where grid folds over itself, or ORCA2 grids
    const double diagonal_NW_SE = diagonal_in_degrees( p_NW, p_SE, sphere_.radius() );
    const double diagonal_SW_NE = diagonal_in_degrees( p_SW, p_NE, sphere_.radius() );
    return std::max( diagonal_NW_SE, diagonal_SW_NE ) > largest_diagonal;
}

bool DetectInvalidElement::diagonal_too_large( const PointLonLat& p_SW, const PointLonLat& p_SE,
                                               const PointLonLat& p_NE, const PointLonLat& p_NW ) const {
    return diagonal_too_large( p_SW, p_SE, p_NE, p_NW, largest_diagonal_ );
}

bool DetectInvalidElement::invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                                            const PointLonLat& p_NW, Statistics& statistics ) const {
    constexpr bool enable_additional_geometry_checks = true;
    constexpr bool enable_southern_edge_aspect_check = true;
    constexpr bool enable_degenerate_element_check = false; // false: allow for degenerate quads... even if not nicely shaped

    constexpr bool enable_quality_metric_check = false;
    constexpr double southern_edge_aspect_margin = 1.1;
    constexpr double minimum_element_quality = 0.1;

    statistics.last_reason = Reason::none;

    PointXYZ xyz_SW{ sphere_.xyz( p_SW ) };
    PointXYZ xyz_SE{ sphere_.xyz( p_SE ) };
    PointXYZ xyz_NE{ sphere_.xyz( p_NE ) };
    PointXYZ xyz_NW{ sphere_.xyz( p_NW ) };

    if (enable_degenerate_element_check) {
        if ( has_zero_edge( xyz_SW, xyz_SE, xyz_NE, xyz_NW, length_tolerance_ ) ) {
            statistics.zero_edge++;
            statistics.invalid_elements++;
            return true;
        }

        if ( has_invalid_jacobian( xyz_SW, xyz_SE, xyz_NE, xyz_NW ) ) {
            statistics.invalid_jacobian++;
            statistics.invalid_elements++;
            return true;
        }
    }
    if ( enable_additional_geometry_checks ) {
        if ( has_self_intersection( xyz_SW, xyz_SE, xyz_NE, xyz_NW ) ) {
            statistics.last_reason = Reason::self_intersecting;
            statistics.self_intersecting++;
            statistics.invalid_elements++;
            return true;
        }

        if ( has_zero_diagonal( xyz_SW, xyz_SE, xyz_NE, xyz_NW, length_tolerance_ ) ) {
            statistics.last_reason = Reason::zero_diagonal;
            statistics.zero_diagonal++;
            statistics.invalid_elements++;
            return true;
        }
    }

    double lat_max = std::max( { p_SW.lat(), p_SE.lat(), p_NE.lat(), p_NW.lat() } );
    if ( invalid_quad_2d( p_SW, p_SE, p_NE, p_NW ) && lat_max < 45. ) {
        statistics.last_reason = Reason::invalid_quad_2d;
        statistics.invalid_quads_2d++;
        statistics.invalid_elements++;
        return true;
    }

    if ( diagonal_too_large( xyz_SW, xyz_SE, xyz_NE, xyz_NW, largest_diagonal_ ) ) {
        statistics.last_reason = Reason::diagonal_too_large;
        statistics.diagonal_too_large++;
        statistics.invalid_elements++;
        return true;
    }

    if ( invalid_quad_3d( xyz_SW, xyz_SE, xyz_NE, xyz_NW ) ) {
        statistics.last_reason = Reason::invalid_quad_3d;
        statistics.invalid_quads_3d++;
        statistics.invalid_elements++;
        return true;
    }

    if ( enable_quality_metric_check ) {
        if ( element_quality( xyz_SW, xyz_SE, xyz_NE, xyz_NW ) < minimum_element_quality ) {
            statistics.last_reason = Reason::poor_quality;
            statistics.poor_quality++;
            statistics.invalid_elements++;
            return true;
        }
    }

    if ( enable_southern_edge_aspect_check && lat_max < -66. ) {
        const double north_south_length =  0.5 * ( chord_length( xyz_SW, xyz_NW ) + chord_length( xyz_SE, xyz_NE ) );
        const double west_east_length   =  0.5 * ( chord_length( xyz_SW, xyz_SE ) + chord_length( xyz_NW, xyz_NE ) );
        if ( north_south_length > southern_edge_aspect_margin * west_east_length ) {
            statistics.last_reason = Reason::southern_edge_aspect;
            statistics.southern_edge_aspect++;
            statistics.invalid_elements++;
            return true;
        }
    }

    if ( eORCA025_ && lat_max < -86.5) {
        statistics.last_reason = Reason::eorca025_southern_cap;
        statistics.invalid_elements++;
        return true;
    }

    if ( ORCA2_ ) {
        if ( lat_max < 60. ) {
            constexpr util::NormaliseLongitude normalized{ -180. };
            double lon_min = normalized( std::min( { p_SW.lon(), p_SE.lon(), p_NE.lon(), p_NW.lon() } ) );
            if ( lon_min > -20. && lon_min < 20. ) {
                double dlon_N = p_NE.lon() - p_NW.lon();
                double dlon_S = p_SE.lon() - p_SW.lon();
                if ( std::max( dlon_N, dlon_S ) > 2. * std::min( dlon_N, dlon_S ) ) {
                    statistics.last_reason = Reason::orca2_longitude_aspect;
                    statistics.invalid_elements++;
                    return true;
                }
            }
        }
        if ( lat_max < 52. && lat_max > 41. ) {
            constexpr util::NormaliseLongitude normalized{ -180. };
            double lon_min = normalized( std::min( { p_SW.lon(), p_SE.lon(), p_NE.lon(), p_NW.lon() } ) );
            if ( lon_min > -2. && lon_min < 4. ) {

                const double north_south_length =
                    0.5 * ( chord_length( xyz_SW, xyz_NW ) + chord_length( xyz_SE, xyz_NE ) );
                const double west_east_length =
                    0.5 * ( chord_length( xyz_SW, xyz_SE ) + chord_length( xyz_NW, xyz_NE ) );
                if ( std::max( north_south_length, west_east_length ) > 2 * std::min(north_south_length, west_east_length ) ) {
                    statistics.last_reason = Reason::western_europe_edge_aspect;
                    statistics.western_europe_edge_aspect++;
                    statistics.invalid_elements++;
                    return true;
                }
            }
        }
    }

    std::array<double,2> diagonals {
        diagonal_in_degrees( xyz_SW, xyz_NE, sphere_.radius() ),
        diagonal_in_degrees( xyz_SE, xyz_NW, sphere_.radius() )
    };
    for ( double diagonal : diagonals ) {
        statistics.minimum_diagonal = std::min( statistics.minimum_diagonal, diagonal );
        statistics.maximum_diagonal = std::max( statistics.maximum_diagonal, diagonal );
        statistics.average_diagonal = ( statistics.average_diagonal * statistics.num_diagonals + diagonal ) / ( statistics.num_diagonals + 1 );
        statistics.num_diagonals++;
    }

    return false;
}

bool DetectInvalidElement::invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE,
                                            const PointLonLat& p_NW ) const {
    Statistics statistics;
    return invalid_element( p_SW, p_SE, p_NE, p_NW, statistics );
}


}  // namespace atlas::orca
