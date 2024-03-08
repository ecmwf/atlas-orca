/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "OrcaGrid.h"

#include "eckit/filesystem/PathName.h"

namespace atlas {

inline const OrcaGrid::grid_t* orca_grid( const Grid::Implementation* grid ) {
    return dynamic_cast<const OrcaGrid::grid_t*>( grid );
}

OrcaGrid::OrcaGrid() : Grid() {}

OrcaGrid::OrcaGrid( const std::string& name ) : Grid{ name }, grid_{ orca_grid( Grid::get() ) } {}

OrcaGrid::OrcaGrid( const Grid& grid ) : Grid{ grid }, grid_{ orca_grid( Grid::get() ) } {}

}  // namespace atlas
