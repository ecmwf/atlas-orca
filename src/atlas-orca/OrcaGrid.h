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

#include <initializer_list>

#include "Orca.h"
#include "atlas/grid/Grid.h"


namespace atlas {

//---------------------------------------------------------------------------------------------------------------------
// Further grid interpretation classes defined in this file

class OrcaGrid;

/*
                                                        Grid
                                                          |
                                    +----------+----------+----------+----------+
                                    |                     |                     |
                             StructuredGrid        UnstructuredGrid          OrcaGrid

*/

//---------------------------------------------------------------------------------------------------------------------

class OrcaGrid : public Grid {
public:
    using grid_t = grid::detail::grid::Orca;

public:
    OrcaGrid();
    OrcaGrid( const std::string& name );
    OrcaGrid( const Grid& );

    bool valid() const { return grid_; }
    operator bool() const { return valid(); }

    idx_t ny() const { return grid_->ny(); }
    idx_t nx() const { return grid_->nx(); }

    using Grid::xy;
    void xy( idx_t i, idx_t j, double xy[] ) const { grid_->xy( i, j, xy ); }

    using Grid::lonlat;
    void lonlat( idx_t i, idx_t j, double lonlat[] ) const { grid_->lonlat( i, j, lonlat ); }

    PointXY xy( idx_t i, idx_t j ) const { return grid_->xy( i, j ); }
    PointLonLat lonlat( idx_t i, idx_t j ) const { return grid_->lonlat( i, j ); }

    double f1( idx_t i, idx_t j ) const { return grid_->lsm( i, j ); }
    double f2( idx_t i, idx_t j ) const { return grid_->core( i, j ); }

    const grid_t* get() const { return grid_; }
    const grid_t* operator->() const { return get(); }

private:
    const grid_t* grid_ = nullptr;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
