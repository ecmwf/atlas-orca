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

#include <bitset>
#include <cstdint>
#include <cstring>

#include "atlas/array/DataType.h"
#include "atlas/runtime/Exception.h"


namespace atlas::orca {

//------------------------------------------------------------------------------------------------------

class Flag {
public:
    static constexpr int GHOST           = 0;
    static constexpr int WATER           = 1;
    static constexpr int INVALID_ELEMENT = 2;

    explicit Flag( std::byte& b ) : byte( b ) { std::memcpy( &bits, &byte, sizeof( std::byte ) ); }
    explicit Flag( const std::byte& b ) : byte( const_cast<std::byte&>( b ) ), read_only{ true } {
        std::memcpy( &bits, &byte, sizeof( std::byte ) );
    }
    void set( int pos ) {
        ATLAS_ASSERT( !read_only );
        bits.set( pos, true );
        std::memcpy( &byte, &bits, sizeof( std::byte ) );
    }
    void unset( int pos ) {
        ATLAS_ASSERT( !read_only );
        bits.set( pos, false );
        std::memcpy( &byte, &bits, sizeof( std::byte ) );
    }
    bool test( int pos ) const { return bits.test( pos ); }

private:
    std::bitset<8> bits;  // 1 byte  ( WARNING, sizeof(std::bitset<8>) != 1; it is compiler dependent, usually 4 or 8 )
    std::byte& byte;
    bool read_only{ false };
};

//------------------------------------------------------------------------------------------------------

}  // namespace atlas::orca
