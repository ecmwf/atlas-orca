/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <limits>

#include "atlas/util/Topology.h"
#include "atlas-orca/meshgenerator/LocalOrcaGrid.h"


#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}
#endif

namespace atlas::orca::meshgenerator {

//----------------------------------------------------------------------------------------------------------------------
LocalOrcaGrid::LocalOrcaGrid(const OrcaGrid& grid, const SurroundingRectangle& rectangle) :
        orca_(grid) {

  // Starting point for local orca grid size is the surrounding rectangle.
  ix_orca_min_ = rectangle.ix_min();
  ix_orca_max_ = rectangle.ix_max();
  iy_orca_min_ = rectangle.iy_min();
  iy_orca_max_ = rectangle.iy_max();

  // Ensure we include the orca halo points if we are at the edge of the orca grid.
  if (rectangle.ix_min() <= 0) {
    ix_orca_min_ = std::min(rectangle.ix_min(), -orca_.haloWest());
  }
  if (rectangle.ix_max() >= orca_.nx() - 1) {
    ix_orca_max_ = std::max(rectangle.ix_max(), orca_.nx() + orca_.haloEast() - 1);
  }
  // no halo at the Southern boundary of orca grid (closed boundary we think?)
  if (rectangle.iy_min() <= 0) {
    iy_orca_min_ = -orca_.haloSouth();
  }

  if (rectangle.iy_max() >= orca_.ny() - 1) {
    iy_orca_max_ = std::max(rectangle.iy_max(), orca_.ny() + orca_.haloNorth() - 1);
  }

  std::cout << " orca_.nx() " << orca_.nx()
            << " orca_.ny() " << orca_.ny() << std::endl;

  std::cout << " ix_orca_min_ " <<  ix_orca_min_
            << " ix_orca_max_ " <<  ix_orca_max_
            << " iy_orca_min_ " <<  iy_orca_min_
            << " iy_orca_max_ " <<  iy_orca_max_ << std::endl;

  // NOTE: +1 because the size of the dimension is one bigger than index of the last element
  nx_orca_ = ix_orca_max_ - ix_orca_min_ + 1;
  ny_orca_ = iy_orca_max_ - iy_orca_min_ + 1;
  size_ = nx_orca_ * ny_orca_;

  std::cout << "nx_orca_ " << nx_orca_ << " ny_orca_ " << ny_orca_ << " size_ " << size_ << std::endl;


  // Partitions and local indices in surrounding rectangle
  parts.resize( size_, -1 );
  halo.resize( size_, 0 );
  is_node.resize( size_, false );
  is_cell.resize(size_, false);
  is_ghost.resize( size_, 1 );
  nb_used_real_nodes_ = 0;
  nb_used_ghost_nodes_ = 0;
  nb_used_real_cells_ = 0;
  nb_used_ghost_cells_ = 0;
  uint16_t nb_used_halo_nodes = 0;
  {
    //atlas_omp_parallel_for( idx_t iy = 0; iy < ny_; iy++ )
    for( size_t iy = 0; iy < ny_orca_; iy++ ) {
      for ( size_t ix = 0; ix < nx_orca_; ix++ ) {
        idx_t ii = index( ix, iy );
        const auto ij_glb_haloed = this->global_ij( ix, iy );
        idx_t ix_reg = ij_glb_haloed.i;
        idx_t iy_reg = ij_glb_haloed.j;
        idx_t reg_ii = -1;
        // Are we in the ORCA halo? If so we use those points instead of wrapping to the
        // other side of the grid.
        if ( this->orca_halo( ix, iy ) ) {
          // Find the information for the closest node to the ORCA halo point
          // that is inside the grid.
          idx_t ix_reg_h = ix_reg;
          idx_t iy_reg_h = iy_reg;
          if ( ix_reg < 0 ) {
            ix_reg_h = 0;
          } else if ( ix_reg >= orca_.nx() ) {
            ix_reg_h = orca_.nx() - 1;
          }
          if ( iy_reg < 0 ) {
            iy_reg_h = 0;
          } else if ( iy_reg >= orca_.ny() ) {
            iy_reg_h = orca_.ny() - 1;
          }
          ix_reg_h -= rectangle.ix_min();
          iy_reg_h -= rectangle.iy_min();
          reg_ii = rectangle.index( ix_reg_h, iy_reg_h );
        } else {
          ix_reg -= rectangle.ix_min();
          iy_reg -= rectangle.iy_min();
          reg_ii = rectangle.index( ix_reg, iy_reg );
        }
        ASSERT(reg_ii >= 0);
        ASSERT(reg_ii < rectangle.parts.size());
        ASSERT(reg_ii < rectangle.halo.size());
        ASSERT(reg_ii < rectangle.is_ghost.size());
        parts.at( ii ) = rectangle.parts.at( reg_ii );
        halo.at( ii ) = rectangle.halo.at( reg_ii );
        is_ghost.at( ii ) = rectangle.is_ghost.at( reg_ii );
      }
    }
  }
  std::cout << "determine number of cells and number of nodes"<< std::endl;

  // determine number of cells and number of nodes
  {
    auto mark_node_used = [&]( int ix, int iy ) {
      idx_t ii = index( ix, iy );
      if ( !is_node.at(ii) ) {
        ++nb_used_nodes_;
        is_node.at(ii) = true;
        if ( is_ghost.at( ii ) ) {
          ++nb_used_ghost_nodes_;
          if ( halo[ii] != 0)
            ++nb_used_halo_nodes;
        } else {
          ++nb_used_real_nodes_;
        }
      }
    };
    auto mark_cell_used = [&]( int ix, int iy ) {
      idx_t ii = index( ix, iy );
      if ( !is_cell.at(ii) ) {
        ++nb_cells_;
        is_cell.at(ii) = true;
        if ( is_ghost.at( ii ) ) {
          ++nb_used_ghost_cells_;
        } else {
          ++nb_used_real_cells_;
        }
      }
    };
    // Loop over all elements to determine which are required
    nb_cells_ = 0;
    nb_used_nodes_ = 0;
    for ( idx_t iy = 0; iy < ny_orca_-1; iy++ ) {
      for ( idx_t ix = 0; ix < nx_orca_-1; ix++ ) {
        idx_t ii = index( ix, iy );
        if ( halo[ii] <= rectangle.halosize() ) {
          mark_cell_used( ix, iy );
          mark_node_used( ix, iy );
          mark_node_used( ix + 1, iy );
          mark_node_used( ix + 1, iy + 1 );
          mark_node_used( ix, iy + 1 );
        }
      }
    }
  }

  // adjust ghost points based on orca halo ghost info
  is_ghost_including_orca_halo.resize( size_, 1 );
  for( size_t iy = 0; iy < ny_orca_; iy++ ) {
    for ( size_t ix = 0; ix < nx_orca_; ix++ ) {
      idx_t ii = index( ix, iy );
      if ( is_node.at( ii ) ) {
        is_ghost_including_orca_halo.at( ii ) = static_cast<bool>(is_ghost.at( ii ));
        const auto ij_glb_haloed = this->orca_haloed_global_grid_ij( ix, iy );
        // The southern boundary does not contain halo points apart from at the
        // east and west limits.
        if ( this->orca_halo( ix, iy ) && 
             ((ij_glb_haloed.j >= 0) || 
              (ij_glb_haloed.i < 0) || 
              (ij_glb_haloed.i >= orca_.nx())) ) {
        // this one should wrap when we have a standard halo, but still
        // preserve the orca halo as though points in the orca are real points
        // for the purposes of the wrapping. But it isn't working quite right
        // The southern boundary does not contain halo points apart from at the
        // east and west limits.
          if ( (ij_glb_haloed.i > orca_.nx() + orca_.haloWest()) ||
               (ij_glb_haloed.j > orca_.ny() + orca_.haloNorth()) ||
               (ij_glb_haloed.i < - orca_.haloEast()) ||
               (ij_glb_haloed.j < - orca_.haloSouth()) ) {
              std::cout << ij_glb_haloed << " out of bounds." << std::endl;
              ASSERT(false);
          }

          is_ghost_including_orca_halo.at( ii ) = static_cast<bool>(is_ghost.at( ii )) || orca_.ghost( ij_glb_haloed.i, ij_glb_haloed.j );
        }
      }
    }
  }

  // setup normalisation objects
  {
    lon00_ = orca_.xy( 0, 0 ).x();
    lon00_normaliser_ = util::NormaliseLongitude(lon00_ - 180. );
  }
}
int LocalOrcaGrid::index( idx_t ix, idx_t iy ) const {
  ATLAS_ASSERT_MSG(static_cast<size_t>(ix) < nx_orca_,
     std::string("ix >= nx_orca_: ") + std::to_string(ix) + " >= " + std::to_string(nx_orca_));
  ATLAS_ASSERT_MSG(static_cast<size_t>(iy) < ny_orca_,
    std::string("iy >= ny_orca_: ") + std::to_string(iy) + " >= " + std::to_string(ny_orca_));
  return iy * nx_orca_ + ix;
}

PointIJ LocalOrcaGrid::global_ij( idx_t ix, idx_t iy ) const {
  return PointIJ(ix_orca_min_ + ix, iy_orca_min_ + iy);
}

const PointXY LocalOrcaGrid::grid_xy( idx_t ix, idx_t iy ) const {
  const auto ij = this->master_global_ij( ix, iy );
  return orca_.xy( ij.i, ij.j );
}

PointXY LocalOrcaGrid::normalised_grid_xy( idx_t ix, idx_t iy ) const {
  double west  = lon00_ - 90.;
  const auto ij = this->master_global_ij( ix, iy );
  if ( (ij.i > orca_.nx() + orca_.haloWest()) ||
       (ij.j > orca_.ny() + orca_.haloNorth()) ||
       (ij.i < - orca_.haloEast()) ||
       (ij.j < - orca_.haloSouth()) ) {
      std::cout << ij << " out of bounds when calling normalised_grid_xy." << std::endl;
      ASSERT(false);
  }
  const PointXY xy = orca_.xy( ij.i, ij.j );
  double lon1 = lon00_normaliser_( orca_.xy( 1, ij.j ).x() );
  if ( lon1 < lon00_ - 10. ) {
      west = lon00_ - 20.;
  }

  if ( ij.i < nx_orca_ / 2 ) {
    auto lon_first_half_normaliser  = util::NormaliseLongitude{west};
    return PointXY( lon_first_half_normaliser( xy.x() ), xy.y() );
  } else {
    auto lon_second_half_normaliser = util::NormaliseLongitude{lon00_ + 90.};
    return PointXY( lon_second_half_normaliser( xy.x() ), xy.y() );
  }
}

// Unique global index of the orca grid excluding orca grid halos (orca halo points are wrapped to their master index).
// Note: right now need to add +1 to this to fill the field with corresponding name.
gidx_t LocalOrcaGrid::master_global_index( idx_t ix, idx_t iy ) const {
  const auto ij = this->global_ij( ix, iy );
  return orca_.periodicIndex( ij.i, ij.j );
}

PointIJ LocalOrcaGrid::master_global_ij( idx_t ix, idx_t iy ) const {
  //const auto ij = this->global_ij( ix, iy );
  //return orca_.periodicIJ( ij.i, ij.j );
  const auto master_idx = this->master_global_index( ix, iy );
  idx_t ix_glb_master, iy_glb_master;
  orca_.index2ij( master_idx, ix_glb_master, iy_glb_master );
  return PointIJ(ix_glb_master, iy_glb_master);
}

PointLonLat LocalOrcaGrid::normalised_grid_master_lonlat( idx_t ix, idx_t iy ) const {
  double west  = lon00_ - 90.;
  const auto ij = this->global_ij( ix, iy );
  if ( (ij.i > orca_.nx() + orca_.haloWest()) ||
       (ij.j > orca_.ny() + orca_.haloNorth()) ||
       (ij.i < - orca_.haloEast()) ||
       (ij.j < - orca_.haloSouth()) ) {
      std::cout << ij << " out of bounds when calling normalised_grid_master_lonlat." << std::endl;
      ASSERT(false);
  }
  //const auto ij = this->orca_haloed_global_grid_ij( ix, iy );
  const auto master_ij = this->master_global_ij( ix, iy );

  const PointLonLat lonlat = orca_.lonlat( master_ij.i, master_ij.j );

  double lon1 = lon00_normaliser_( orca_.xy( 1, ij.j ).x() );
  if ( lon1 < lon00_ - 10. ) {
      west = lon00_ - 20.;
  }

  if ( ij.i < nx_orca_ / 2 ) {
    auto lon_first_half_normaliser  = util::NormaliseLongitude{west};
    return PointLonLat( lon_first_half_normaliser( lonlat.lon() ), lonlat.lat() );
  } else {
    auto lon_second_half_normaliser = util::NormaliseLongitude{lon00_ + 90.};
    return PointLonLat( lon_second_half_normaliser( lonlat.lon() ), lonlat.lat() );
  }
}

// Unique global index of the orca grid including orca grid halos.  This needs
// to wrap points back into the orca grid, but is subtly different from
// OrcaGrid.periodicIndex as it will only wrap points that are outside of the
// orca grid halos
PointIJ LocalOrcaGrid::orca_haloed_global_grid_ij( idx_t ix, idx_t iy ) const {
  // global grid properties
  const auto iy_glb_min = -orca_.haloSouth();
  const auto ix_glb_min = -orca_.haloWest();
  const auto nx_orca_glb = orca_.nx() + orca_.haloEast() + orca_.haloWest();
  const auto ny_orca_glb = orca_.ny() + orca_.haloSouth() + orca_.haloNorth();
  const idx_t glbarray_offset  = -( nx_orca_glb * iy_glb_min ) - ix_glb_min;
  const idx_t glbarray_jstride = nx_orca_glb;

  auto ij = this->global_ij( ix, iy );

  ATLAS_ASSERT_MSG( ij.j >= -orca_.haloSouth(),
                    std::string("Cannot extend below the southern ORCA halo: j = ") +
                    std::to_string(ij.j) + std::string(", orca_.haloSouth() = ") + 
                    std::to_string(orca_.haloSouth()) );

  // wrap points outside of orca_grid halo back into the orca grid.
  if ( (ij.i >= ix_glb_min + nx_orca_glb) || (ij.j >= iy_glb_min + ny_orca_glb)
    || (ij.i < ix_glb_min) ) {
    gidx_t p_idx = orca_.periodicIndex(ij.i, ij.j);
    idx_t i, j;
    orca_.index2ij(p_idx, i, j);
    ij.i = i;
    ij.j = j;
    //ij = orca_.periodicIJ(ij.i, ij.j);
  }

  ATLAS_ASSERT_MSG( ij.i > ix_glb_min,
                    std::to_string(ij.i) + std::string(" < ") + std::to_string(ix_glb_min) );
  ATLAS_ASSERT_MSG( ij.j > iy_glb_min,
                    std::to_string(ij.j) + std::string(" < ") + std::to_string(iy_glb_min) );
  ATLAS_ASSERT_MSG( ij.i < ix_glb_min + nx_orca_glb,
                    std::to_string(ij.i) + std::string(" > ") + std::to_string(ix_glb_min + nx_orca_glb) );
  ATLAS_ASSERT_MSG( ij.j < iy_glb_min + ny_orca_glb,
                    std::to_string(ij.j) + std::string(" > ") + std::to_string(iy_glb_min + ny_orca_glb) );
  return ij;
}

idx_t LocalOrcaGrid::orca_haloed_global_grid_index( idx_t ix, idx_t iy ) const {
  PointIJ ij = this->orca_haloed_global_grid_ij( ix, iy );
  const auto iy_glb_min = -orca_.haloSouth();
  const auto ix_glb_min = -orca_.haloWest();
  const auto nx_orca_glb = orca_.nx() + orca_.haloEast() + orca_.haloWest();
  const auto ny_orca_glb = orca_.ny() + orca_.haloSouth() + orca_.haloNorth();
  const idx_t glbarray_offset  = -( nx_orca_glb * iy_glb_min ) - ix_glb_min;
  const idx_t glbarray_jstride = nx_orca_glb;
  return glbarray_offset + ij.j * glbarray_jstride + ij.i;
}

void LocalOrcaGrid::flags( idx_t ix, idx_t iy, util::detail::BitflagsView<int>& flag_view ) const {
  flag_view.reset();

  const auto ij_glb = this->orca_haloed_global_grid_ij( ix, iy );

  if ( (ij_glb.i > orca_.nx() + orca_.haloWest()) ||
       (ij_glb.j > orca_.ny() + orca_.haloNorth()) ||
       (ij_glb.i < - orca_.haloEast()) ||
       (ij_glb.j < - orca_.haloSouth()) ) {
      std::cout << ij_glb << " out of bounds when calling orca_.flags." << std::endl;
      ASSERT(false);
  }

  if ( this->is_ghost[this->index(ix, iy)] ) {
    flag_view.set( util::Topology::GHOST );
    if( this->orca_haloed_global_grid_index( ix, iy ) !=
        this->master_global_index( ix, iy ) ) {
      const auto normalised_xy = this->normalised_grid_xy( ix, iy );
      if ( ij_glb.i >= orca_.nx() - orca_.haloWest() ) {
        flag_view.set( util::Topology::PERIODIC );
      }
      else if ( ij_glb.i < orca_.haloEast() - 1 ) {
        flag_view.set( util::Topology::PERIODIC );
      }
      if ( ij_glb.j >= orca_.ny() - orca_.haloNorth() - 1 ) {
          flag_view.set( util::Topology::PERIODIC );
          if ( normalised_xy.x() > lon00_ + 90. ) {
              flag_view.set( util::Topology::EAST );
          }
          else {
              flag_view.set( util::Topology::WEST );
          }
      }

      if ( flag_view.check( util::Topology::PERIODIC ) ) {
        // It can still happen that nodes were flagged as periodic wrongly
        // e.g. where the grid folds into itself
        auto lonlat_master = this->normalised_grid_master_lonlat(ix, iy);
        if( std::abs(lonlat_master.lon() - normalised_xy.x()) < 1.e-12 ) {
            flag_view.unset( util::Topology::PERIODIC );
            // if (( std::abs(lonlat_master.lat() - normalised_xy.y()) < 1.e-12 ) &&
            //     ( ij_glb.j >= orca_.ny() - orca_.haloNorth() - 1 ))
            //     grid_fold_inodes.push_back(inode);
        }
      }
    }
  }

  flag_view.set( orca_.land( ij_glb.i, ij_glb.j ) ? util::Topology::LAND : util::Topology::WATER );

  if ( ij_glb.i <= 0 ) {
    flag_view.set( util::Topology::BC | util::Topology::WEST );
  }
  else if ( ij_glb.j >= orca_.nx() ) {
    flag_view.set( util::Topology::BC | util::Topology::EAST );
  }
}
bool LocalOrcaGrid::water( idx_t ix, idx_t iy ) const {

  const auto ij_glb = this->orca_haloed_global_grid_ij( ix, iy );

  const auto ij_glb = this->global_ij( ix, iy );
  if ( (ij_glb.i > orca_.nx() + orca_.haloWest()) ||
       (ij_glb.j > orca_.ny() + orca_.haloNorth()) ||
       (ij_glb.i < - orca_.haloEast()) ||
       (ij_glb.j < - orca_.haloSouth()) ) {
      std::cout << ij_glb << " out of bounds when calling orca_.water." << std::endl;
      ASSERT(false);
  }

  return orca_.water( ij_glb.i, ij_glb.j );
}
bool LocalOrcaGrid::orca_halo( idx_t ix, idx_t iy ) const {
  const auto ij_glb = this->orca_haloed_global_grid_ij( ix, iy );
  if ( ((ij_glb.i < 0) && (ij_glb.i >= -orca_.haloWest())) ||
       ((ij_glb.i >= orca_.nx()) && (ij_glb.i < (orca_.nx() + orca_.haloEast()))) ||
       ((ij_glb.j < 0) && (ij_glb.j >= -orca_.haloSouth())) ||
       ((ij_glb.j >= orca_.ny()) && (ij_glb.j < (orca_.ny() + orca_.haloNorth()))) ) {
    return true;
  }
  return false;
}
}  // namespace atlas::orca::meshgenerator
