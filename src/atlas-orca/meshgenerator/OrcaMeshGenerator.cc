/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "OrcaMeshGenerator.h"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "eckit/utils/Hash.h"

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/StructuredMeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Topology.h"

#include "atlas-orca/meshgenerator/SurroundingRectangle.h"
#include "atlas-orca/meshgenerator/LocalOrcaGrid.h"


namespace atlas::orca::meshgenerator {

// ORCA2 interesting indices
static const std::vector<std::pair<int, int>> glb_ij_pairs{
    {0, 148},
    {1, 148},
    {0, 147},
    {89, 148},
    {90, 148},
    {91, 148},
    {92, 148},
    {93, 148},
    //89-90, 147
    {89, 147},
    {90, 147},
    //92-101, 147
    {92, 147},
    {93, 147},
    {94, 147},
    {95, 147},
    {96, 147},
    {97, 147},
    {98, 147},
    {99, 147},
    {100, 147},
    {101, 147},
    //89-90, 146
    {89, 146},
    {90, 146},
    //180-181, 148
    {180, 148},
    {181, 148},
    {181, 146},
};

namespace {

StructuredGrid equivalent_regular_grid( const OrcaGrid& orca ) {
    ATLAS_ASSERT( orca );
    // Mimic hole in South pole, and numbering from South to North. patch determines if endpoint is at North Pole
    StructuredGrid::YSpace yspace{grid::LinearSpacing{{-80., 90.}, orca.ny(), true}};  //not patch.at( orca.name() )}};
    // Periodic xspace
    StructuredGrid::XSpace xspace{grid::LinearSpacing{{0., 360.}, orca.nx(), false}};

    return StructuredGrid{xspace, yspace};
}
}  // namespace

struct Nodes {
    array::ArrayView<idx_t, 2> ij;
    array::ArrayView<double, 2> xy;
    array::ArrayView<double, 2> lonlat;
    array::ArrayView<gidx_t, 1> glb_idx;
    array::IndexView<idx_t, 1> remote_idx;
    array::ArrayView<int, 1> part;
    array::ArrayView<int, 1> ghost;
    array::ArrayView<int, 1> halo;
    array::ArrayView<int, 1> node_flags;
    array::ArrayView<int, 1> water;
    array::ArrayView<gidx_t, 1> master_glb_idx;

    util::detail::BitflagsView<int> flags( idx_t i ) { return util::Topology::view( node_flags( i ) ); }

    explicit Nodes( Mesh& mesh ) :
        ij{ array::make_view<idx_t, 2>( mesh.nodes().add(
            Field( "ij", array::make_datatype<idx_t>(), array::make_shape( mesh.nodes().size(), 2 ) ) ) ) },
        xy{ array::make_view<double, 2>( mesh.nodes().xy() ) },
        lonlat{ array::make_view<double, 2>( mesh.nodes().lonlat() ) },
        glb_idx{ array::make_view<gidx_t, 1>( mesh.nodes().global_index() ) },
        remote_idx{ array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() ) },
        part{ array::make_view<int, 1>( mesh.nodes().partition() ) },
        ghost{ array::make_view<int, 1>( mesh.nodes().ghost() ) },
        halo{ array::make_view<int, 1>( mesh.nodes().halo() ) },
        node_flags{ array::make_view<int, 1>( mesh.nodes().flags() ) },
        water{ array::make_view<int, 1>( mesh.nodes().add(
            Field( "water", array::make_datatype<int>(), array::make_shape( mesh.nodes().size() ) ) ) ) },
        master_glb_idx{ array::make_view<gidx_t, 1>( mesh.nodes().add( Field(
            "master_global_index", array::make_datatype<gidx_t>(), array::make_shape( mesh.nodes().size() ) ) ) ) } {}
};

struct Cells {
    array::ArrayView<int, 1> part;
    array::ArrayView<int, 1> halo;
    array::ArrayView<gidx_t, 1> glb_idx;
    array::ArrayView<int, 1> flags_view;
    mesh::HybridElements::Connectivity& node_connectivity;
    util::detail::BitflagsView<int> flags( idx_t i ) { return util::Topology::view( flags_view( i ) ); }
    explicit Cells( Mesh& mesh ) :
        part{ array::make_view<int, 1>( mesh.cells().partition() ) },
        halo{ array::make_view<int, 1>( mesh.cells().halo() ) },
        glb_idx{ array::make_view<gidx_t, 1>( mesh.cells().global_index() ) },
        flags_view{ array::make_view<int, 1>( mesh.cells().flags() ) },
        node_connectivity( mesh.cells().node_connectivity() ) {}
};

void OrcaMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    ATLAS_TRACE( "OrcaMeshGenerator::generate" );
    using Topology = util::Topology;

    OrcaGrid orca_grid{grid};
    ATLAS_ASSERT( orca_grid );
    ATLAS_ASSERT( !mesh.generated() );

    // global (all processor) configuration information about ORCA grid for the ORCA mesh under construction
    SurroundingRectangle::Configuration SR_cfg;
    SR_cfg.mypart = mypart_;
    SR_cfg.nparts = nparts_;
    SR_cfg.halosize = halosize_;
    SR_cfg.nx_glb = orca_grid.nx();
    SR_cfg.ny_glb = orca_grid.ny();

    SurroundingRectangle SR(distribution, SR_cfg);
    LocalOrcaGrid local_orca(orca_grid, SR);

    // global orca grid dimensions and index limits
    auto ny_orca_halo = orca_grid.ny() + orca_grid.haloNorth() + orca_grid.haloSouth();
    auto nx_orca_halo = orca_grid.nx() + orca_grid.haloEast() + orca_grid.haloWest();
    auto iy_glb_min = -orca_grid.haloSouth();
    auto iy_glb_max = orca_grid.ny() + orca_grid.haloNorth() - 1;
    auto ix_glb_min = -orca_grid.haloWest();
    auto ix_glb_max = orca_grid.nx() + orca_grid.haloEast() - 1;

    // clone some grid properties
    setGrid( mesh, grid, distribution );

    const bool serial_distribution = (SR_cfg.nparts == 1 || distribution.type() == "serial");
    if ( serial_distribution && (halosize_ > 0) ) {
      throw_NotImplemented("Halo size should be 0 if using a serial distribution or a single partition", Here());
    }

    //---------------------------------------------------

    if ( serial_distribution ) {
        ATLAS_ASSERT_MSG(ix_glb_max == local_orca.ix_max(),
          std::string("Size of the surrounding rectangle x-space doesn't match up with orca-grid x-space: ")
          + std::to_string(ix_glb_max) + " != " + std::to_string(local_orca.ix_max()) );
        ATLAS_ASSERT_MSG(iy_glb_max == local_orca.iy_max(),
          std::string("Size of the surrounding rectangle y-space doesn't match up with orca-grid y-space: ")
          + std::to_string(iy_glb_max) + " != " + std::to_string(local_orca.iy_max()) );
        ATLAS_ASSERT_MSG(ix_glb_min == local_orca.ix_min(),
          std::string("Size of the surrounding rectangle x-space doesn't match up with orca-grid x-space: ")
          + std::to_string(ix_glb_min) + " != " + std::to_string(local_orca.ix_min()) );
        ATLAS_ASSERT_MSG(iy_glb_min == local_orca.iy_min(),
          std::string("Size of the surrounding rectangle y-space doesn't match up with orca-grid y-space: ")
          + std::to_string(iy_glb_min) + " != " + std::to_string(local_orca.iy_min()) );
        ATLAS_ASSERT_MSG(nx_orca_halo * ny_orca_halo == local_orca.nb_used_nodes(),
          std::string("Number of used nodes doesn't match total size of orca grid plus halo: ")
          + std::to_string(nx_orca_halo * ny_orca_halo) + " != " + std::to_string(local_orca.nb_used_nodes()) );
    }
    ATLAS_ASSERT_MSG(local_orca.nb_used_real_nodes() + local_orca.nb_used_ghost_nodes() == local_orca.nb_used_nodes(),
      std::string("Number of nodes in mesh does not equal the total size of the surrounding rectangle")
      + std::to_string(local_orca.nb_used_real_nodes() + local_orca.nb_used_ghost_nodes()) + " != " + std::to_string(local_orca.nb_used_nodes()) );

    // define nodes and associated properties
    mesh.nodes().resize(local_orca.nb_used_nodes());
    Nodes nodes( mesh );

    // define cells and associated properties
#if ATLAS_TEMPORARY_ELEMENTTYPES
    // DEPRECATED
    mesh.cells().add( new mesh::temporary::Quadrilateral(), local_orca.nb_cells() );
#else
    // Use this since atlas 0.35.0
    mesh.cells().add( mesh::ElementType::create("Quadrilateral"), local_orca.nb_cells() );
#endif

    Cells cells( mesh );

    int inode_nonghost = 0;
    int inode_ghost    = 0;

    int ix_pivot = SR_cfg.nx_glb / 2;
    bool patch   = not orca_grid.ghost( ix_pivot + 1, SR_cfg.ny_glb - 1 );

    std::vector<idx_t> node_index( local_orca.nx()*local_orca.ny(), -1 );

    {
        ATLAS_TRACE( "nodes" );

        // loop over nodes and set properties
        inode_nonghost = 0;
        inode_ghost    = local_orca.nb_used_real_nodes();  // orca ghost nodes don't count as ghost nodes for this stuff?... something to do with the orca grid?

        ATLAS_TRACE_SCOPE( "indexing" )
        for ( idx_t iy = 0; iy < local_orca.ny(); iy++ ) {
            for ( idx_t ix = 0; ix < local_orca.nx(); ix++ ) {
                idx_t ii = local_orca.index( ix, iy );
                // node properties
                if ( local_orca.is_node[ii] ) {
                    // set node counter
                    ATLAS_ASSERT_MSG( ii < node_index.size(),
                        std::to_string(ii) + " >= " + std::to_string(node_index.size()));
                    if ( local_orca.is_ghost[ii] ) {
                        node_index.at(ii) = inode_ghost++;
                    } else {
                        node_index.at(ii) = inode_nonghost++;
                    }
                    ATLAS_ASSERT_MSG( node_index[ii] < local_orca.nb_used_nodes(),
                        std::string("node_index[") + std::to_string(ii) + std::string("] ") + std::to_string(node_index[ii])
                        + " >= " + std::to_string(local_orca.nb_used_nodes()));
                }
            }
        }

        ATLAS_TRACE_SCOPE( "filling" )
        for( idx_t iy = 0; iy < local_orca.ny(); iy++ ) {
            for ( idx_t ix = 0; ix < local_orca.nx(); ix++ ) {
                idx_t ii = local_orca.index( ix, iy );
                // node properties
                if ( local_orca.is_node[ii] ) {
                    idx_t inode = node_index[ii];
                    ASSERT(ii < local_orca.is_ghost.size());
                    ASSERT(ii < local_orca.nx() * local_orca.ny());
                    ASSERT(inode < local_orca.nb_used_nodes());

                    // ghost nodes
                    nodes.ghost( inode ) = local_orca.is_ghost_including_orca_halo[ii];

                    // global index
                    nodes.glb_idx( inode ) = local_orca.orca_haloed_global_grid_index( ix, iy ) + 1;

                    // grid ij coordinates
                    {
                      const auto ij_glb = local_orca.orca_haloed_global_grid_ij( ix, iy );
                      nodes.ij( inode, XX ) = ij_glb.i;
                      nodes.ij( inode, YY ) = ij_glb.j;
                    }

                    const auto normalised_xy = local_orca.normalised_grid_xy( ix, iy );

                    // grid xy coordinates
                    {
                      nodes.xy( inode, LON ) = normalised_xy.x();
                      nodes.xy( inode, LAT ) = normalised_xy.y();
                    }

                    // geographic coordinates (normalised)
                    {
                      nodes.lonlat( inode, LON ) = normalised_xy.x();
                      nodes.lonlat( inode, LAT ) = normalised_xy.y();
                    }

                    // part and remote_idx
                    nodes.part( inode )           = local_orca.parts[ii];
                    nodes.remote_idx( inode )     = inode;
                    nodes.master_glb_idx( inode ) = nodes.glb_idx( inode );

                    // flags
                    auto flags = nodes.flags( inode );

                    if ( nodes.ghost( inode ) ) {
                      gidx_t master_idx             = local_orca.master_global_index( ix, iy );
                      nodes.master_glb_idx( inode ) = master_idx + 1;
                      if ( nparts_ == 1 ) {
                        nodes.part( inode ) = 0;
                      } else {
                        PointIJ master_ij = local_orca.master_global_ij( ix, iy );
                        auto clamp = []( idx_t value, idx_t lower, idx_t upper ) {
                          // in C++17 this is std::clamp
                          return std::max( lower, std::min( value, upper ) );
                        };
                        idx_t clamped_i = clamp( master_ij.i, 0, SR_cfg.nx_glb - 1 );
                        idx_t clamped_j = clamp( master_ij.j, 0, SR_cfg.ny_glb - 1 );
                        nodes.part( inode ) = SR.global_partition( clamped_i, clamped_j );
                      }
                      nodes.remote_idx( inode ) = serial_distribution ?
                          static_cast<int>( master_idx ) : -1;
                    }

                    if ( serial_distribution ) {
                      ATLAS_ASSERT(nodes.remote_idx( inode ) > -1,
                                   std::string("remote index is -1:  " + std::to_string(nodes.remote_idx( inode ))));
                    }

                    local_orca.flags( ix, iy, flags );

                    nodes.water( inode ) = local_orca.water( ix, iy );
                    nodes.halo( inode ) = local_orca.halo[ii];
                }
            }
        }
    }
    std::vector<idx_t> cell_index( local_orca.nx()*local_orca.ny() );
    // loop over nodes and define cells, putting non-ghost cells first
    {
        ATLAS_TRACE( "elements" );
        idx_t jcell = 0;
        idx_t jcell_ghost = local_orca.nb_used_real_cells();
        ATLAS_TRACE_SCOPE( "indexing" );
        for ( idx_t iy = 0; iy < local_orca.ny() - 1; iy++ ) {      // don't loop into ghost/periodicity row
            for ( idx_t ix = 0; ix < local_orca.nx() - 1; ix++ ) {  // don't loop into ghost/periodicity column
                idx_t ii = local_orca.index( ix, iy );
                std::stringstream assert_msg;
                assert_msg << ii << " > " << cell_index.size() << std::endl;
                ATLAS_ASSERT(ii <  cell_index.size(), assert_msg.str());
                if ( local_orca.is_ghost[ii] ) {
                    cell_index[ii] = jcell_ghost++;
                } else {
                    cell_index[ii] = jcell++;
                }
            }
        }

        ATLAS_TRACE_SCOPE( "filling" )
        //atlas_omp_parallel_for( idx_t iy = 0; iy < local_orca.ny() - 1; iy++ ) {  // don't loop into ghost/periodicity row
        for( idx_t iy = 0; iy < local_orca.ny() - 1; iy++ ) {  // don't loop into ghost/periodicity row
            for ( idx_t ix = 0; ix < local_orca.nx() - 1; ix++ ) {                // don't loop into ghost/periodicity column
                idx_t ii   = local_orca.index( ix, iy );
                int ix_glb = local_orca.ix_min() + ix;
                int iy_glb = local_orca.iy_min() + iy;
                if ( local_orca.is_cell[ii] ) {
                    idx_t jcell = cell_index[ii];

                    // define cell corners (local indices)
                    std::array<idx_t, 4> quad_nodes{};
                    quad_nodes[0] = node_index[local_orca.index( ix, iy )];          // lower left
                    quad_nodes[1] = node_index[local_orca.index( ix + 1, iy )];      // lower right
                    quad_nodes[2] = node_index[local_orca.index( ix + 1, iy + 1 )];  // upper right
                    quad_nodes[3] = node_index[local_orca.index( ix, iy + 1 )];      // upper left

                    cells.flags( jcell ).reset();

                    cells.node_connectivity.set( jcell, quad_nodes.data() );
                    cells.part( jcell )    = nodes.part( quad_nodes[0] );
                    cells.glb_idx( jcell ) = ( iy_glb - iy_glb_min ) * nx_orca_halo + ( ix_glb - ix_glb_min ) + 1;

                    if ( iy_glb >= SR_cfg.ny_glb - 1 ) {
                        cells.flags( jcell ).set( Topology::GHOST );
                        if ( patch && ix_glb < ix_pivot ) {                     // case of eg ORCA1_T
                            cells.part( jcell ) = nodes.part( quad_nodes[0] );  // lower left
                            cells.flags( jcell ).unset( Topology::GHOST );
                        }
                        else {                                                  // case of eg ORCA2_T
                            cells.part( jcell ) = nodes.part( quad_nodes[2] );  // upper right
                        }
                    }

                    bool elem_contains_water_point = [&] {
                        for ( idx_t inode : quad_nodes ) {
                            if ( nodes.flags( inode ).check( Topology::WATER ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    bool elem_contains_land_point = [&] {
                        for ( idx_t inode : quad_nodes ) {
                            if ( nodes.flags( inode ).check( Topology::LAND ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    cells.halo( jcell ) = [&] {
                        int h = 0;
                        for ( idx_t inode : quad_nodes ) {
                            h = std::max( h, nodes.halo( inode ) );
                        }
                        if ( iy_glb < 0 ) {
                            h = 0;
                            if ( ix_glb < 0 ) {
                                h = -ix_glb;
                            }
                            else if ( ix_glb >= SR_cfg.nx_glb ) {
                                h = ix_glb - ( SR_cfg.nx_glb - 1 );
                            }
                        }
                        return h;
                    }();

                    if ( elem_contains_water_point ) {
                        cells.flags( jcell ).set( Topology::WATER );
                    }
                    if ( elem_contains_land_point ) {
                        cells.flags( jcell ).set( Topology::LAND );
                    }
                    const auto ij_glb = local_orca.orca_haloed_global_grid_ij( ix, iy );
                    if ( orca_grid.invalidElement( ij_glb.i, ij_glb.j ) ) {
                        cells.flags( jcell ).set( Topology::INVALID );
                    }
                }
            }
        }
    }
    ATLAS_DEBUG_VAR( serial_distribution );
    if ( serial_distribution ) {
        // Bypass for "BuildParallelFields"
        mesh.nodes().metadata().set( "parallel", true );

        // Bypass for "BuildPeriodicBoundaries"
        mesh.metadata().set( "periodic", true );
    }
    else {
        ATLAS_DEBUG( "build_remote_index" );
        build_remote_index( mesh );
    }

    // Degenerate points in the ORCA mesh mean that the standard BuildHalo
    // methods for updating halo sizes will not work.
    mesh.metadata().set("halo_locked", true);
    mesh.metadata().set("halo", halosize_);
    mesh.nodes().metadata().set<size_t>( "NbRealPts", local_orca.nb_used_nodes() );
    mesh.nodes().metadata().set<size_t>( "NbVirtualPts", 0 );
}

using Unique2Node = std::map<gidx_t, idx_t>;
void OrcaMeshGenerator::build_remote_index(Mesh& mesh) const {
    ATLAS_TRACE();

    mesh::Nodes& nodes = mesh.nodes();

    bool parallel = false;
    bool periodic = false;
    nodes.metadata().get( "parallel", parallel );
    mesh.metadata().get( "periodic", periodic );
    if ( parallel || periodic ) {
        ATLAS_DEBUG( "build_remote_index: already parallel, return" );
        return;
    }

    auto mpi_size = mpi::size();
    auto mypart   = mpi::rank();
    int nb_nodes  = nodes.size();

    // get the indices and partition data
    auto master_glb_idx = array::make_view<gidx_t, 1>( nodes.field( "master_global_index" ) );
    auto glb_idx        = array::make_view<gidx_t, 1>( nodes.global_index() );
    auto ridx           = array::make_indexview<idx_t, 1>( nodes.remote_index() );
    auto part           = array::make_view<int, 1>( nodes.partition() );
    auto ghost          = array::make_view<int, 1>( nodes.ghost() );
    auto ij             = array::make_view<idx_t, 2>( nodes.field("ij"));

    // find the nodes I want to request the data for
    std::vector<std::vector<gidx_t>> send_uid( mpi_size );
    std::vector<std::vector<int>> req_lidx( mpi_size );

    Unique2Node global2local;
    for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
        gidx_t uid = master_glb_idx( jnode );
        if ( ghost( jnode ) ) {
            send_uid[part( jnode )].push_back( uid );
            req_lidx[part( jnode )].push_back( jnode );
            ridx( jnode ) = -1;
        }
        else {
            ridx( jnode ) = jnode;
        }
        if ( ghost( jnode ) == 0 ) {
            bool inserted = global2local.insert( std::make_pair( uid, jnode ) ).second;
            ATLAS_ASSERT( inserted, std::string( "index already inserted " ) + std::to_string( uid ) + ", " +
                                        std::to_string( jnode ) + " at jnode " + std::to_string( global2local[uid] ) );
        }
    }

    std::vector<std::vector<gidx_t>> recv_uid( mpi_size );

    // Request data from those indices
    mpi::comm().allToAll( send_uid, recv_uid );

    // Find and populate send vector with indices to send
    std::vector<std::vector<int>> send_ridx( mpi_size );
    std::vector<std::vector<int>> send_gidx( mpi_size );
    std::vector<std::vector<int>> send_part( mpi_size );
    for ( idx_t p = 0; p < mpi_size; ++p ) {
        for ( idx_t i = 0; i < recv_uid[p].size(); ++i ) {
            idx_t found_idx = -1;
            gidx_t uid      = recv_uid[p][i];
            if ( auto found = global2local.find( uid ); found != global2local.end() ) {
                found_idx = found->second;
            }

            ATLAS_ASSERT( found_idx != -1, "master global index not found: " + std::to_string( uid ) + " for sent element index: " + std::to_string( i ) );
            send_ridx[p].push_back( ridx( found_idx ) );
            send_gidx[p].push_back( static_cast<int>( glb_idx( found_idx ) ) );
            send_part[p].push_back( part( found_idx ) );
        }
    }

    std::vector<std::vector<int>> recv_ridx( mpi_size );
    std::vector<std::vector<int>> recv_gidx( mpi_size );
    std::vector<std::vector<int>> recv_part( mpi_size );

    mpi::comm().allToAll( send_ridx, recv_ridx );
    mpi::comm().allToAll( send_gidx, recv_gidx );
    mpi::comm().allToAll( send_part, recv_part );

    // Fill out missing remote indices
    for ( idx_t p = 0; p < mpi_size; ++p ) {
        for ( idx_t i = 0; i < recv_ridx[p].size(); ++i ) {
            ridx( req_lidx[p][i] )    = recv_ridx[p][i];
            glb_idx( req_lidx[p][i] ) = recv_gidx[p][i];
            part( req_lidx[p][i] )    = recv_part[p][i];
        }
    }

    // sanity check
    for ( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
        ATLAS_ASSERT( ridx( jnode ) >= 0, "ridx not filled with part " + std::to_string( part( jnode ) ) + " at " +
                                              std::to_string( jnode ) );
    }

    mesh.metadata().set( "periodic", true );
    nodes.metadata().set( "parallel", true );
}

OrcaMeshGenerator::OrcaMeshGenerator( const eckit::Parametrisation& config ) {
    config.get( "partition", mypart_ = mpi::rank() );
    config.get( "partitions", nparts_ = mpi::size() );
    config.get( "halo", halosize_);
    if ( halosize_ < 0 ) {
      throw_NotImplemented("Halo size must be >= 0", Here());
    }
}

void OrcaMeshGenerator::generate( const Grid& grid, const grid::Partitioner& partitioner, Mesh& mesh ) const {
    std::unordered_set<std::string> valid_distributions = {"serial", "checkerboard", "equal_regions", "equal_area"};
    ATLAS_ASSERT(valid_distributions.find(partitioner.type()) != valid_distributions.end(),
                 partitioner.type() + " is not an implemented distribution type. "
                 + "Valid types are 'serial', 'checkerboard' or 'equal_regions', 'equal_area'");
    if (partitioner.type() == "serial" && halosize_ >= 1)
      throw_NotImplemented("halo size must be zero for 'serial' distribution type ORCA grids", Here());
    auto regular_grid = equivalent_regular_grid( grid );
    auto distribution = grid::Distribution( regular_grid, partitioner );
    generate( grid, distribution, mesh );
}

void OrcaMeshGenerator::generate( const Grid& grid, Mesh& mesh ) const {
    generate( grid, grid::Partitioner( grid.partitioner() ), mesh );
}

void OrcaMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "OrcaMeshGenerator" );
}

namespace {
atlas::meshgenerator::MeshGeneratorBuilder<OrcaMeshGenerator> __OrcaMeshGenerator( "orca" );
}

}  // namespace atlas::orca::meshgenerator
