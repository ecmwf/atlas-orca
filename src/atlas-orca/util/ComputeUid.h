/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace atlas {
namespace orca {

//------------------------------------------------------------------------------------------------------

std::string compute_uid( const std::string& arrangement, const double lon[], const double lat[], int nx_halo,
                         int ny_halo, const std::int32_t halo[] );

std::string compute_uid( const std::string& arrangement, const std::vector<double>& lon, const std::vector<double>& lat,
                         int nx_halo, int ny_halo, const std::array<std::int32_t, 4>& halo );

//------------------------------------------------------------------------------------------------------

}  // namespace orca
}  // namespace atlas
