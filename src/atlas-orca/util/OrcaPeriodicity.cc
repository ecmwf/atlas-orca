/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "OrcaPeriodicity.h"

#include <array>
#include <cstdint>

#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/OrcaData.h"
#include "atlas-orca/util/PointIJ.h"
#include "atlas/library/config.h"
namespace atlas {
namespace orca {

OrcaPeriodicity::OrcaPeriodicity( const OrcaData& orca ) : OrcaPeriodicity( orca.dimensions, orca.halo, orca.pivot ) {}

OrcaPeriodicity::OrcaPeriodicity( const std::array<std::int32_t, 2>& dimensions,
                                  const std::array<std::int32_t, 4>& halo, const std::array<double, 2>& pivot ) :
    nx_{dimensions[0] - halo_[HALO_WEST] - halo_[HALO_EAST]},
    ny_{dimensions[0] - halo_[HALO_WEST] - halo_[HALO_EAST]},
    halo_{halo},
    pivot_{pivot} {
    ibegin_ = halo_[HALO_WEST];
    iend_   = ibegin_ + nx_;
}

//------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------

}  // namespace orca
}  // namespace atlas
