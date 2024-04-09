/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "NetCDFReader.h"

#if ATLAS_ORCA_HAVE_NETCDF
#include <netcdf>
#endif

#include <algorithm>
#include <sstream>

#include "eckit/types/Fraction.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/vector.h"

#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"


namespace atlas::orca {

NetCDFReader::NetCDFReader( const util::Config& config ) {
    config.get( "arrangement", arrangement_ );
    std::string P = arrangement_;
    if ( P != "T" && P != "W" && P != "F" && P != "U" && P != "V" ) {
        ATLAS_THROW_EXCEPTION( "arrangement expected to be any of {F,T,U,V,W}. Received: " << P );
    }
}

void NetCDFReader::read( const std::string& uri, OrcaData& data ) {
#if ATLAS_ORCA_HAVE_NETCDF == 0
    ATLAS_THROW_EXCEPTION( "atlas-orca was not compiled with NetCDF support" );
#else
    OrcaDataFile file{ uri };

    netCDF::NcFile ncdata( file, netCDF::NcFile::read );
    auto read_dimensions = [&]( size_t& ni, size_t& nj ) {
        auto get_dim = [&]( const std::string& name ) {
            auto dim = ncdata.getDim( name );
            if ( dim.isNull() ) {
                ATLAS_THROW_EXCEPTION( "Dimension '" << name << "' is not present in NetCDF file" );
            }
            return dim.getSize();
        };
        ni = get_dim( "x" );
        nj = get_dim( "y" );
    };

    auto read_variable = [&]( const std::string& name, std::vector<double>& values ) {
        auto var = ncdata.getVar( name );
        if ( var.isNull() ) {
            ATLAS_THROW_EXCEPTION( "Variable '" << name << "' is not present in NetCDF file" );
        }
        size_t ni = 0;
        size_t nj = 0;
        read_dimensions( ni, nj );
        auto dims = var.getDims();
        values.resize( ni * nj );
        if ( dims.size() == 3 ) {
            ATLAS_ASSERT( dims[1].getSize() == nj );
            ATLAS_ASSERT( dims[2].getSize() == ni );
            var.getVar( { 0, 0, 0 }, { 1, nj, ni }, values.data() );
        }
        else if ( dims.size() == 4 ) {
            ATLAS_ASSERT( dims[2].getSize() == nj );
            ATLAS_ASSERT( dims[3].getSize() == ni );
            var.getVar( { 0, 0, 0, 0 }, { 1, 1, nj, ni }, values.data() );
        }
        else {
            std::stringstream errmsg;
            errmsg << "Unexpected dimensions for variable '" << name << "': Found [";
            for ( size_t i = 0; i < dims.size(); ++i ) {
                errmsg << dims[i].getSize();
                errmsg << ( i < dims.size() - 1 ? "," : "]" );
            }
            ATLAS_THROW_EXCEPTION( errmsg.str() );
        }
    };

    size_t ni     = 0;
    size_t nj     = 0;
    std::string P = arrangement_;
    ATLAS_ASSERT( P.size() == 1 );
    std::string p = P;
    p[0]          = ( P == "W" ) ? 't' : static_cast<char>( std::tolower( P[0] ) );
    std::vector<double> mask;

    read_dimensions( ni, nj );
    read_variable( "glam" + p, data.lon );
    read_variable( "gphi" + p, data.lat );
    read_variable( p + "mask", mask );

    data.dimensions = { static_cast<int>( ni ), static_cast<int>( nj ) };
    data.flags.resize( ni * nj );
    for ( size_t n = 0; n < data.flags.size(); ++n ) {
        Flag flag{ data.flags[n] };
        if ( mask[n] != 0. ) {
            flag.set( Flag::WATER );
        }
    }

    data.halo[HALO_EAST] = 1;
    data.halo[HALO_WEST] = 1;
    eckit::Fraction resolution{ 360, static_cast<int>( ni ) - data.halo[HALO_EAST] - data.halo[HALO_WEST] };
    if ( resolution != 1. ) {
        // T-pivot
        data.pivot = { static_cast<double>( ni / 2 ), static_cast<double>( nj - 2 ) };
        if ( P == "T" || P == "W" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 1;
        }
        else if ( P == "F" ) {
            data.halo[HALO_NORTH] = 2;
            data.halo[HALO_SOUTH] = 0;
            data.pivot[0] -= 0.5;
            data.pivot[1] -= 0.5;
        }
        else if ( P == "U" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 1;
            data.pivot[0] -= 0.5;
        }
        else if ( P == "V" ) {
            data.halo[HALO_NORTH] = 2;
            data.halo[HALO_SOUTH] = 1;
            data.pivot[1] -= 0.5;
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }
    else {
        // F-pivot
        data.pivot = { static_cast<double>( ni / 2 - 1 ), static_cast<double>( nj - 2 ) };
        if ( P == "T" || P == "W" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 1;
            data.pivot[0] += 0.5;
            data.pivot[1] += 0.5;
        }
        else if ( P == "F" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 0;
        }
        else if ( P == "U" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 1;
            data.pivot[1] += 0.5;
        }
        else if ( P == "V" ) {
            data.halo[HALO_NORTH] = 1;
            data.halo[HALO_SOUTH] = 1;
            data.pivot[0] += 0.5;
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }

    data.setGhost();
    data.makeHaloConsistent();
#endif
}

}  // namespace atlas::orca
