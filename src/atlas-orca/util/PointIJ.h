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

#include <iostream>

#include "atlas/library/config.h"


namespace atlas::orca {

//------------------------------------------------------------------------------------------------------

struct PointIJ {
    PointIJ() = default;
    PointIJ( idx_t _i, idx_t _j ) : i( _i ), j( _j ) {}
    idx_t i;
    idx_t j;
    friend std::ostream& operator<<( std::ostream& out, const PointIJ& p ) {
        out << "{" << p.i << "," << p.j << "}";
        return out;
    }
    bool operator==( const PointIJ& other ) const { return other.i == i && other.j == j; }
    bool operator!=( const PointIJ& other ) const { return other.i != i || other.j != j; }
    bool operator<( const PointIJ& other ) const {
        if ( j < other.j ) {
            return true;
        }
        if ( j == other.j && i < other.i ) {
            return true;
        }
        return false;
    }
    bool operator>( const PointIJ& other ) const {
        if ( j > other.j ) {
            return true;
        }
        if ( j == other.j && i > other.i ) {
            return true;
        }
        return false;
    }
};

//------------------------------------------------------------------------------------------------------

}  // namespace atlas::orca
