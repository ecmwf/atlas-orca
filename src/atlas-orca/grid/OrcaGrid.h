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

#include <initializer_list>

#include "Orca.h"
#include "atlas/grid/Grid.h"
#include "atlas-orca/util/PointIJ.h"


namespace atlas {

class OrcaGrid : public Grid {
public:
    using grid_t = grid::detail::grid::Orca;

public:
    OrcaGrid();
    OrcaGrid( const std::string& name );
    OrcaGrid( const Grid& );

    bool valid() const { return grid_ != nullptr; }
    explicit operator bool() const { return valid(); }

    idx_t ny() const { return grid_->ny(); }
    idx_t nx() const { return grid_->nx(); }

    using Grid::xy;
    void xy( idx_t i, idx_t j, double xy[] ) const { grid_->xy( i, j, xy ); }

    using Grid::lonlat;
    void lonlat( idx_t i, idx_t j, double lonlat[] ) const { grid_->lonlat( i, j, lonlat ); }

    PointXY xy( idx_t i, idx_t j ) const { return grid_->xy( i, j ); }
    PointLonLat lonlat( idx_t i, idx_t j ) const { return grid_->lonlat( i, j ); }

    bool water( idx_t i, idx_t j ) const { return grid_->water( i, j ); }
    bool land( idx_t i, idx_t j ) const { return grid_->land( i, j ); }
    bool ghost( idx_t i, idx_t j ) const { return grid_->ghost( i, j ); }
    bool invalidElement( idx_t i, idx_t j ) const { return grid_->invalidElement( i, j ); }

    int haloWest() const { return grid_->haloWest(); }
    int haloEast() const { return grid_->haloEast(); }
    int haloNorth() const { return grid_->haloNorth(); }
    int haloSouth() const { return grid_->haloSouth(); }

    gidx_t periodicIndex( idx_t i, idx_t j ) const { return grid_->periodicIndex( i, j ); }
    orca::PointIJ periodicIJ( idx_t i, idx_t j ) const { return grid_->periodicIJ( i, j ); }

    void index2ij( gidx_t gidx, idx_t& i, idx_t& j ) const { grid_->index2ij( gidx, i, j ); }

    const grid_t* get() const { return grid_; }
    const grid_t* operator->() const { return get(); }

private:
    const grid_t* grid_ = nullptr;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
