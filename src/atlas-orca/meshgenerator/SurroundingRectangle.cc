/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SurroundingRectangle.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <fstream>
#include <iomanip>

#include "atlas/array/Array.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Geometry.h"


namespace atlas::orca::meshgenerator {

namespace {
int wrap( idx_t value, idx_t lower, idx_t upper ) {
  // wrap around coordinate system when out of bounds
  const idx_t width = upper - lower;
  if (value < lower) {
    return wrap(value + width, lower, upper);
  }
  if (value >= upper) {
    return wrap(value - width, lower, upper);
  }
  return value;
}
}  // namespace

PointIJ SurroundingRectangle::global_periodic_ij(idx_t ix_glb, idx_t iy_glb) const {
  // wrap around coordinate system when out of bounds on the global rectangle

  // use global ij on rectangle
  const idx_t width_x = cfg_.nx_glb;
  const idx_t ix_glb_max = cfg_.nx_glb - 1;
  const idx_t iy_glb_max = cfg_.ny_glb - 1;
  idx_t ix_glb_p = ix_glb;
  idx_t iy_glb_p = iy_glb;

  // j index north/south boundaries
  if (iy_glb_p > iy_glb_max) {
    ix_glb_p = ix_glb_p + width_x/2;
    iy_glb_p = 2*iy_glb_max - iy_glb_p;
  }
  if (iy_glb_p < 0) {
    ix_glb_p = ix_glb_p + width_x/2;
    iy_glb_p = -iy_glb_p;
  }

  // i index periodic east/west boundaries
  if (ix_glb_p < 0) {
    ix_glb_p = wrap(ix_glb_p + width_x, 0, cfg_.nx_glb);
  }
  if (ix_glb_p > ix_glb_max) {
    ix_glb_p = wrap(ix_glb_p - width_x, 0, cfg_.nx_glb);
  }

  // convert to local ij on rectangle
  return PointIJ(ix_glb_p, iy_glb_p);
}

int SurroundingRectangle::index( int ix, int iy ) const {
  ATLAS_ASSERT_MSG(ix < nx_, std::string("ix >= nx_: ") + std::to_string(ix) + " >= " + std::to_string(nx_));
  ATLAS_ASSERT_MSG(iy < ny_, std::string("iy >= ny_: ") + std::to_string(iy) + " >= " + std::to_string(ny_));
  return iy * nx_ + ix;
}

int SurroundingRectangle::partition( idx_t i, idx_t j ) const {
  PointIJ ij = this->global_periodic_ij(ix_min_ + i, iy_min_ + j);
  ATLAS_ASSERT_MSG(ij.i < cfg_.nx_glb, std::string("ix >= cfg_.nx_glb: ") + std::to_string(ij.i) + " >= " + std::to_string(cfg_.nx_glb) + " ix_min_ " + std::to_string(ix_min_));
  ATLAS_ASSERT_MSG(ij.j < cfg_.ny_glb, std::string("iy >= cfg_.ny_glb: ") + std::to_string(ij.j) + " >= " + std::to_string(cfg_.ny_glb) + " iy_min_ " + std::to_string(iy_min_));
  return distribution_.partition( ij.j * cfg_.nx_glb + ij.i );
}

int SurroundingRectangle::global_partition( idx_t ix_glb, idx_t iy_glb ) const {
  PointIJ ij = this->global_periodic_ij(ix_glb, iy_glb);
  auto ix_glb_p = ij.i;
  auto iy_glb_p = ij.j;
  ATLAS_ASSERT_MSG(ix_glb_p < cfg_.nx_glb, std::string("ix >= cfg_.nx_glb: ") + std::to_string(ix_glb_p) + " >= " + std::to_string(cfg_.nx_glb));
  ATLAS_ASSERT_MSG(iy_glb_p < cfg_.ny_glb, std::string("iy >= cfg_.ny_glb: ") + std::to_string(iy_glb_p) + " >= " + std::to_string(cfg_.ny_glb));
  return distribution_.partition( iy_glb_p * cfg_.nx_glb + ix_glb_p );
}

SurroundingRectangle::SurroundingRectangle(
    const grid::Distribution& distribution,
    const Configuration& cfg )
  : distribution_(distribution), cfg_(cfg) {
  ATLAS_TRACE();
  cfg_.check_consistency();

  // determine rectangle (ix_min_:ix_max_) x (iy_min_:iy_max_) surrounding the nodes on this processor
  ix_min_         = cfg_.nx_glb;
  ix_max_         = 0;
  iy_min_         = cfg_.ny_glb;
  iy_max_         = 0;
  nb_real_nodes_owned_by_rectangle = 0;

  {
    ATLAS_TRACE( "find rectangle bounds" );
    atlas_omp_parallel {
      int ix_min_TP = ix_min_;
      int ix_max_TP = ix_max_;
      int iy_min_TP = iy_min_;
      int iy_max_TP = iy_max_;
      int nb_real_nodes_owned_by_rectangle_TP = 0;
      atlas_omp_for( idx_t iy_glb = 0; iy_glb < cfg_.ny_glb; iy_glb++ ) {
        for ( idx_t ix_glb = 0; ix_glb < cfg_.nx_glb; ix_glb++ ) {
          ATLAS_ASSERT_MSG(ix_glb < cfg_.nx_glb, std::string("ix_glb >= cfg_.nx_glb: ") + std::to_string(ix_glb) + " >= " + std::to_string(cfg_.nx_glb));
          ATLAS_ASSERT_MSG(iy_glb < cfg_.ny_glb, std::string("iy_glb >= cfg_.ny_glb: ") + std::to_string(iy_glb) + " >= " + std::to_string(cfg_.ny_glb));
          int p = global_partition( ix_glb, iy_glb );
          if ( p == cfg_.mypart ) {
            ix_min_TP = std::min<idx_t>( ix_min_TP, ix_glb );
            ix_max_TP = std::max<idx_t>( ix_max_TP, ix_glb );
            iy_min_TP = std::min<idx_t>( iy_min_TP, iy_glb );
            iy_max_TP = std::max<idx_t>( iy_max_TP, iy_glb );
            nb_real_nodes_owned_by_rectangle_TP++;
          }
        }
      }
      atlas_omp_critical {
        nb_real_nodes_owned_by_rectangle += nb_real_nodes_owned_by_rectangle_TP;
        ix_min_ = std::min<int>( ix_min_TP, ix_min_);
        ix_max_ = std::max<int>( ix_max_TP, ix_max_);
        iy_min_ = std::min<int>( iy_min_TP, iy_min_);
        iy_max_ = std::max<int>( iy_max_TP, iy_max_);
      }
    }
  }

  // add the halo.
  halosize_ = cfg_.halosize;
  ix_min_ -= halosize_;
  ix_max_ += halosize_;
  iy_min_ -= halosize_;
  iy_max_ += halosize_;

  // +1 to surround the ghost nodes used to complete the cells
  ix_max_ += 1;
  iy_max_ += 1;

  // dimensions of the surrounding rectangle (+1 because the size of the dimension is one bigger than the index of the last element)
  nx_ = ix_max_ - ix_min_ + 1;
  ny_ = iy_max_ - iy_min_ + 1;

  // upper estimate for number of nodes
  uint64_t size = ny_ * nx_;

  // partitions and local indices in SR
  parts.resize( size, -1 );
  halo.resize( size, 0 );
  is_ghost.resize( size, true );
  // vectors marking nodes that are necessary for this proc's cells

  {
    ATLAS_TRACE( "partition, is_ghost, halo" );
    //atlas_omp_parallel_for( idx_t iy = 0; iy < ny_; iy++ )
    for( idx_t iy = 0; iy < ny_; iy++ ) {
      for ( idx_t ix = 0; ix < nx_; ix++ ) {
        idx_t ii = index( ix, iy );
        parts.at( ii ) = partition( ix, iy );
        
        PointIJ pij = global_periodic_ij( ix_min_ + ix, iy_min_ + iy );
        bool periodic_point = ( (pij.i != (ix_min_ + ix) ) || (pij.j != (iy_min_ + iy)) );
        bool halo_found = false;
        int halo_dist = halosize_;
        if ((halosize_ > 0) && ((parts.at( ii ) != cfg_.mypart) || periodic_point) ) {
          // search the surrounding halosize index square for a node on my
          // partition to determine the halo distance
          for (idx_t searchsize = 1; std::max(nx_, ny_); ++searchsize) {
            for (idx_t dhy = -searchsize; dhy <= searchsize; ++dhy) {
              for (idx_t dhx = -searchsize; dhx <= searchsize; ++dhx) {
                if ( ((std::abs(dhx) != searchsize) && (std::abs(dhy) != searchsize)) ||
                     (ix + dhx < 0) || (ix + dhx >= nx_) ||
                     (iy + dhy < 0) || (iy + dhy >= ny_) ) {
                  continue;
                }
                if (partition(ix + dhx, iy + dhy) == cfg_.mypart) {
                  // find the minimum distance from this halo node to
                  // a node on the partition
                  auto dist = std::max(std::abs(dhx), std::abs(dhy));
                  halo_dist = std::min(dist, halo_dist);
                  halo_found = true;
                }
              }
            }
            if ( halo_found ) {
                break;
            }
          }
          ATLAS_ASSERT_MSG( halo_found, std::string("Halo distance not found at point ") +
                                        std::to_string(ix) + std::string(", ") + 
                                        std::to_string(iy) ); 
          halo.at( ii ) = halo_dist;
        }
        is_ghost.at( ii ) = ( (parts.at( ii ) != cfg_.mypart) || periodic_point );
      }
    }
  }
}

}  // namespace atlas::orca::meshgenerator 
