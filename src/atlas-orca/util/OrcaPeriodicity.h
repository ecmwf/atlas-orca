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

#include <cstdint>
#include <array>

#include "atlas/library/config.h"
#include "atlas-orca/util/PointIJ.h"

namespace atlas {
namespace orca {

//------------------------------------------------------------------------------------------------------

class OrcaData;

class OrcaPeriodicity {
private:
    std::array<double,2> pivot_;
    std::array<std::int32_t,4> halo_;
    idx_t nx_;
    idx_t ny_;
    idx_t ibegin_;
    idx_t iend_;

public:
    OrcaPeriodicity( const OrcaData& orca );

    OrcaPeriodicity( const std::array<std::int32_t,2>& dimensions, const std::array<std::int32_t,4>& halo, const std::array<double,2>& pivot );

    PointIJ operator()(PointIJ p ) const {
        return compute(p.i,p.j);
    }

    PointIJ operator()( idx_t i, idx_t j ) const {
        return compute(i,j);
    }

    PointIJ compute( idx_t i, idx_t j ) const {
        PointIJ master{i,j};
        if( i < ibegin_ ) {
            master.i += nx_;
        }
        if( i >= iend_ ) {
            master.i -= nx_;
        }
        if( double(master.j) > pivot_[1] || (double(master.j) == pivot_[1] && double(master.i) > pivot_[0]) ) {
            master.i = idx_t(2.*pivot_[0]) - master.i;
            master.j = idx_t(2.*pivot_[1]) - master.j;
            if( master.i < ibegin_ ) {
                master.i += nx_;
            }
            if( master.i >= iend_ ) {
                master.i -= nx_;
            }
        }
        return master;
    }
};

}  // namespace orca
}  // namespace atlas
