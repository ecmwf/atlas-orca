#pragma once

#include <algorithm>
#include <numeric>
#include <tuple>
#include <utility>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Topology.h"
#include "atlas/util/Point.h"
#include "atlas/util/Config.h"
#include "atlas/interpolation/element/Quad3D.h"

namespace atlas {
namespace orca {
namespace meshgenerator {

class DetectInvalidElement {
public:

    DetectInvalidElement( const util::Config& config ) {
        config.get("ORCA2",orca2_);
        config.get("diagonal",largest_diatonal_);
    }

    bool invalid_quad_2d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE, const PointLonLat& p_NW ) {
        double dlat_W = p_NW.lat() - p_SW.lat();
        double dlat_E = p_NE.lat() - p_SE.lat();
        double dlon_N = p_NE.lon() - p_NW.lon();
        double dlon_S = p_SE.lon() - p_SW.lon();
        return (dlat_W < -1.e-10 || dlat_E < -1.e-10 || dlon_N < -1.e-10 || dlon_S < -1.e-10 );
    }

    bool invalid_quad_3d( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE, const PointLonLat& p_NW ) {
        PointXYZ xyz_SW { sphere_.xyz(p_SW) };
        PointXYZ xyz_SE { sphere_.xyz(p_SE) };
        PointXYZ xyz_NE { sphere_.xyz(p_NE) };
        PointXYZ xyz_NW { sphere_.xyz(p_NW) };
        return not atlas::interpolation::element::Quad3D{xyz_SW,xyz_SE,xyz_NE,xyz_NW}.validate();
    }

    bool diagonal_too_large(  const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE, const PointLonLat& p_NW, double largest_diagonal) {
        if( largest_diagonal == 0. ) {
            return false;
        }
        // example with extended grids where grid folds over itself, or ORCA2 grids
        double diagonal_treshold_2 = largest_diagonal*largest_diagonal;
        double d2_NW_SE = PointLonLat::distance2(p_NW,p_SE);
        double d2_SW_NE = PointLonLat::distance2(p_SW,p_NE);
        return std::max(d2_NW_SE,d2_SW_NE) > diagonal_treshold_2;
    }

    bool diagonal_too_large(  const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE, const PointLonLat& p_NW) {
        return diagonal_too_large(p_SW,p_SE,p_NE,p_NW,largest_diatonal_);
    }

    bool invalid_element( const PointLonLat& p_SW, const PointLonLat& p_SE, const PointLonLat& p_NE, const PointLonLat& p_NW ) {
        double lat_max = std::max({p_SW.lat(),p_SE.lat(),p_NE.lat(),p_NW.lat()});

        if( invalid_quad_2d(p_SW,p_SE,p_NE,p_NW) && lat_max < 45. ) {
            return true;
        }
        if( lat_max < 60. && diagonal_too_large(p_SW,p_SE,p_NE,p_NW) ) {
            return true;
        }
        if( invalid_quad_3d(p_SW,p_SE,p_NE,p_NW) ) {
            return true;
        }
        if( orca2_ ) {
            if( lat_max < 60. ) {
                constexpr util::NormaliseLongitude normalized{-180.};
                double lon_min = normalized( std::min({p_SW.lon(),p_SE.lon(),p_NE.lon(),p_NW.lon()}) );
                if( lon_min > -20. && lon_min < 20. ) {
                    double dlat_W = p_NW.lat() - p_SW.lat();
                    double dlat_E = p_NE.lat() - p_SE.lat();
                    double dlon_N = p_NE.lon() - p_NW.lon();
                    double dlon_S = p_SE.lon() - p_SW.lon();
                    if( std::max(dlon_N,dlon_S) > 2.* std::min(dlon_N,dlon_S) ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }


private:
    geometry::Earth sphere_;
    double largest_diatonal_{0};
    bool orca2_{false};
};

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
