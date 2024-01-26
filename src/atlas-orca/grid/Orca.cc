/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Orca.h"


#include <algorithm>
#include <cctype>

#include "eckit/utils/Hash.h"

#include "atlas/domain/Domain.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

#include "atlas-orca/Library.h"
#include "atlas-orca/util/AtlasIOReader.h"
#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"
#include "atlas-orca/util/OrcaData.h"
#include "atlas-orca/util/OrcaDataFile.h"
#include "atlas-orca/util/OrcaPeriodicity.h"
#include "atlas-orca/util/PointIJ.h"


namespace atlas::grid::detail::grid {

static bool validate_uid() {
    static bool ATLAS_ORCA_VALIDATE_UID = bool(
        eckit::LibResource<bool, atlas::orca::Library>( "atlas-orca-validate-uid;$ATLAS_ORCA_VALIDATE_UID", false ) );
    return ATLAS_ORCA_VALIDATE_UID;
}

static bool compute_uid() {
    static bool ATLAS_ORCA_COMPUTE_UID = bool(
        eckit::LibResource<bool, atlas::orca::Library>( "atlas-orca-compute-uid;$ATLAS_ORCA_COMPUTE_UID", false ) );
    return ATLAS_ORCA_COMPUTE_UID;
}

static std::string spec_name( const Grid::Spec& spec, const std::string& def = "" ) {
    auto n = spec.getString( "orca_name", "" );
    auto a = spec.getString( "orca_arrangement", "" );
    return n.empty() || a.empty() ? def : ( n + "_" + a );
}

static std::string spec_uid( const Grid::Spec& spec, const std::string& def = "" ) {
    return spec.getString( "uid", def );
}

namespace {

template <typename T>
size_t memory( const T& container ) {
    return sizeof( typename T::value_type ) * container.size();
}
size_t memory( const std::vector<bool>& container ) {
    return container.size() / 8;
}

}  // namespace

std::string Orca::static_type() {
    return "ORCA";
}

std::string Orca::name() const {
    return name_;
}

std::string Orca::type() const {
    return static_type();
}

gidx_t Orca::periodicIndex( idx_t i, idx_t j ) const {
    auto p = periodicity_->compute( i + imin_, j + jmin_ );
    return index( p.i - imin_, p.j - jmin_ );
}

orca::PointIJ Orca::periodicIJ( idx_t i, idx_t j ) const {
    return periodicity_->compute( i + imin_, j + jmin_ );
}

Orca::Orca( const Config& config ) :
    Orca( spec_name( config, spec_uid( config, config.getString( "name", "" ) ) ), config ) {}

Orca::Orca( const std::string& name, const Config& config ) :
    name_( spec_name( config, spec_uid( config, name ) ) ),
    imin_( 0 ),
    jmin_( 0 ),
    istride_( 0 ),
    jstride_( 0 ),
    spec_( config ) {
    auto trace = atlas::Trace( Here(), "Orca(" + name_ + ")" );

    orca::OrcaData data;

    {
        // Read data
        orca::OrcaDataFile file{ spec_.getString( "data" ) };
        orca::AtlasIOReader{}.read( file, data );
    }

    halo_north_ = data.halo[orca::HALO_NORTH];
    halo_west_  = data.halo[orca::HALO_WEST];
    halo_east_  = data.halo[orca::HALO_EAST];
    halo_south_ = data.halo[orca::HALO_SOUTH];
    nx_halo_    = data.dimensions[0];
    ny_halo_    = data.dimensions[1];
    nx_         = nx_halo_ - halo_west_ - halo_east_;
    ny_         = ny_halo_ - halo_south_ - halo_north_;
    imin_       = halo_west_;
    jmin_       = halo_south_;
    istride_    = 1;
    jstride_    = nx_halo_;
    points_.resize( nx_halo_ * ny_halo_ );
    water_.resize( nx_halo_ * ny_halo_ );
    ghost_.resize( nx_halo_ * ny_halo_ );
    invalid_element_.resize( nx_halo_ * ny_halo_ );

    size_t nn{ 0 };
    for ( idx_t j = 0; j < ny_halo_; ++j ) {
        for ( idx_t i = 0; i < nx_halo_; ++i, ++nn ) {
            gidx_t n = j * jstride_ + i;
            std::bitset<8> bits;
            std::memcpy( &bits, &data.flags[n], sizeof( std::byte ) );
            ghost_[nn]           = bits.test( orca::Flag::GHOST );
            water_[nn]           = bits.test( orca::Flag::WATER );
            invalid_element_[nn] = bits.test( orca::Flag::INVALID_ELEMENT );
            points_[nn].assign( data.lon[n], data.lat[n] );
        }
    }

    if ( validate_uid() ) {
        Log::debug() << "Validating uid of ORCA grid" << std::endl;
        std::string computed_uid = data.computeUid( spec_ );
        std::string uid          = spec_.getString( "uid" );
        if ( uid != computed_uid ) {
            ATLAS_THROW_EXCEPTION( "ORCA grid encoded uid does not validate with computed uid: " << computed_uid );
        }
    }
    else if ( compute_uid() ) {
        Log::debug() << "Computing uid of ORCA grid" << std::endl;
        spec_.set( "uid", data.computeUid( spec_ ) );
    }

    domain_      = GlobalDomain();
    periodicity_ = std::make_unique<orca::OrcaPeriodicity>( data );
}

void Orca::hash( eckit::Hash& h ) const {
    h.add( uid() );
}

std::string Orca::uid() const {
    return spec_uid( spec_ );
}

idx_t Orca::size() const {
    return nx_halo_ * ny_halo_;
}

RectangularLonLatDomain Orca::lonlatBoundingBox() const {
    return domain_;
}

Grid::Spec Orca::spec() const {
    return spec_;
}

void Orca::print( std::ostream& os ) const {
    os << "Orca(" << name_ << ")";
}

size_t Orca::footprint() const {
    return memory( points_ ) + memory( water_ ) + memory( ghost_ );
}

Grid::Config Orca::meshgenerator() const {
    return { "type", "orca" };
}

Grid::Config Orca::partitioner() const {
    return Config( "type", "checkerboard" )( "regular", true );
}

}  // namespace atlas::grid::detail::grid
