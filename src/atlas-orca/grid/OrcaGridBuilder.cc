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
#include <iomanip>
#include <iostream>

#include "atlas-orca/Library.h"

#include "atlas/grid/SpecRegistry.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/library.h"
#include "atlas/util/Config.h"


namespace atlas::grid::detail::grid {


static class OrcaGridBuilder : public GridBuilder {
    using Implementation = atlas::Grid::Implementation;
    using Config         = Grid::Config;

public:
    OrcaGridBuilder() :
        GridBuilder( Orca::static_type(), { "^e?ORCA[0-9]+_[FTUVW]$" }, { "[e]ORCA<deg>_{F,T,U,V,W}" } ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 30 ) << "[e]ORCA<deg>_{F,T,U,V,W}"
           << "ORCA Tripolar grid. Possible increasing resolutions <deg>: 2,1,025,12";
    }

    const Implementation* create( const std::string& name_or_uid, const Config& /* config */ ) const override {
        auto sane_id( name_or_uid );
        std::transform( sane_id.begin(), sane_id.end(), sane_id.begin(), ::tolower );

        if ( SpecRegistry::has( sane_id ) ) {
            return create( SpecRegistry::get( sane_id ) );
        }

        auto sane_name( name_or_uid );
        std::transform( sane_name.begin(), sane_name.end(), sane_name.begin(), ::toupper );
        if ( sane_name.front() == 'E' ) {
            sane_name.front() = 'e';
        }

        if ( SpecRegistry::has( sane_name ) ) {
            return create( SpecRegistry::get( sane_name ) );
        }

        return nullptr;
    }

    const Implementation* create( const Config& config ) const override {
        std::string type;
        config.get( "type", type );
        if ( type != "ORCA" ) {
            return nullptr;
        }
        return new Orca( config );
    }

    void force_link() {}

} orca_;

}  // namespace atlas::grid::detail::grid
