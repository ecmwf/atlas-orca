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
#include <utility>

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

#include "atlas-orca/grid/OrcaGrid.h"


namespace atlas {
namespace orca {
namespace meshgenerator {


struct Configuration {
    int nparts;
    int mypart;
};

struct SurroundingRectangle {
    std::vector<int> parts;
    std::vector<int> is_ghost;
    // WARNING vector<bool> is not thread-safe, don't assign in parallel region
    std::vector<bool> is_node;
    int size;
    int nx, ny;
    int ix_min, ix_max;
    int iy_min, iy_max;
    int nb_nodes;
    int nb_cells;
    int nb_nodes_owned;

    int index( int i, int j ) const { return j * nx + i; }
    int index_sp( int i ) const { return nx * ny + i; }

    SurroundingRectangle( const Grid& grid, const grid::Distribution& distribution, const Configuration& cfg ) {
        ATLAS_TRACE();
        OrcaGrid orca{grid};
        int mypart     = cfg.mypart;
        int nx_glb     = orca.nx();
        int ny_glb     = orca.ny();
        int ny_halo    = ny_glb + orca.haloNorth();
        int iy_glb_max = ny_glb + orca.haloNorth() - 1;
        int iy_glb_min = -orca.haloSouth();
        int ix_glb_max = nx_glb + orca.haloEast() - 1;
        int ix_glb_min = -orca.haloWest();

        auto partition = [&]( idx_t i, idx_t j ) -> int {
            auto clamp = []( idx_t value, idx_t lower, idx_t upper ) {
                // in C++17 this is std::clamp
                return std::max( lower, std::min( value, upper ) );
            };
            i = clamp( i, 0, nx_glb - 1 );
            j = clamp( j, 0, ny_glb - 1 );
            return distribution.partition( j * nx_glb + i );
        };

        // determine rectangle (ix_min:ix_max) x (iy_min:iy_max) surrounding the nodes on this processor
        ix_min         = nx_glb + 1;
        ix_max         = 0;
        iy_min         = ny_glb + 1;
        iy_max         = 0;
        nb_nodes_owned = 0;

        {
            ATLAS_TRACE( "bounds" );
            atlas_omp_parallel {
                int ix_min_TP = ix_min;
                int ix_max_TP = ix_max;
                int iy_min_TP = iy_min;
                int iy_max_TP = iy_max;
                int nb_nodes_owned_TP = 0;
                atlas_omp_for( idx_t iy = iy_glb_min; iy <= iy_glb_max; iy++ ) {
                    for ( idx_t ix = ix_glb_min; ix <= ix_glb_max; ix++ ) {
                        int p = partition( ix, iy );
                        if ( p == mypart ) {
                            ix_min_TP = std::min<idx_t>( ix_min_TP, ix );
                            ix_max_TP = std::max<idx_t>( ix_max_TP, ix );
                            iy_min_TP = std::min<idx_t>( iy_min_TP, iy );
                            iy_max_TP = std::max<idx_t>( iy_max_TP, iy );
                            nb_nodes_owned_TP++;
                        }
                    }
                }
                atlas_omp_critical {
                    nb_nodes_owned += nb_nodes_owned_TP;
                    ix_min = std::min<int>( ix_min_TP, ix_min);
                    ix_max = std::max<int>( ix_max_TP, ix_max);
                    iy_min = std::min<int>( iy_min_TP, iy_min);
                    iy_max = std::max<int>( iy_max_TP, iy_max);
                }
            }
        }

        // add one row/column for ghost nodes (which include periodicity points)
        ix_max = std::min( ix_glb_max, ix_max + 1 );
        iy_max = std::min( iy_glb_max, iy_max + 1 );

        // dimensions of surrounding rectangle (SR)
        nx = ix_max - ix_min + 1;
        ny = iy_max - iy_min + 1;

        // upper estimate for number of nodes
        size = nx * ny;
        ;

        // partitions and local indices in SR
        parts.resize( size, -1 );
        is_ghost.resize( size, true );
        // vectors marking nodes that are necessary for this proc's cells
        is_node.resize( size, false );

        {  // Compute : SR.part, SR.is_ghost
            ATLAS_TRACE( "part,is_ghost" );
            atlas_omp_parallel_for( idx_t iy = 0; iy < ny; iy++ ) {
                idx_t iy_glb = iy_min + iy;  // global y-index
                for ( idx_t ix = 0; ix < nx; ix++ ) {
                    idx_t ii       = index( ix, iy );
                    idx_t ix_glb   = ix_min + ix;  // global x-index
                    parts.at( ii ) = partition( ix_glb, iy_glb );

                    if ( iy_glb < ny_halo ) {
                        is_ghost.at( ii ) = ( parts.at( ii ) != mypart );
                    }
                }
            }
        }
        // determine number of cells and number of nodes
        {  // Compute SR.is_node
            ATLAS_TRACE( "is_node" );
            auto mark_node_used = [&]( int ix, int iy ) {
                idx_t ii = index( ix, iy );
                if ( !is_node[ii] ) {
                    ++nb_nodes;
                    is_node[ii] = true;
                }
            };
            // Loop over all elements
            nb_cells = 0;
            nb_nodes = 0;
            for ( idx_t iy = 0; iy < ny - 1; iy++ ) {      // don't loop into ghost/periodicity row
                for ( idx_t ix = 0; ix < nx - 1; ix++ ) {  // don't loop into ghost/periodicity column
                    if ( !is_ghost[index( ix, iy )] ) {
                        ++nb_cells;
                        mark_node_used( ix, iy );
                        mark_node_used( ix + 1, iy );
                        mark_node_used( ix + 1, iy + 1 );
                        mark_node_used( ix, iy + 1 );
                    }
                }
            }
        }
    }
};

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
    //array::ArrayView<int, 1> core;
    array::ArrayView<gidx_t, 1> master_glb_idx;

    util::detail::BitflagsView<int> flags( idx_t i ) { return util::Topology::view( node_flags( i ) ); }

    Nodes( Mesh& mesh ) :
        ij{array::make_view<idx_t, 2>( mesh.nodes().add(
            Field( "ij", array::make_datatype<idx_t>(), array::make_shape( mesh.nodes().size(), 2 ) ) ) )},
        xy{array::make_view<double, 2>( mesh.nodes().xy() )},
        lonlat{array::make_view<double, 2>( mesh.nodes().lonlat() )},
        glb_idx{array::make_view<gidx_t, 1>( mesh.nodes().global_index() )},
        remote_idx{array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() )},
        part{array::make_view<int, 1>( mesh.nodes().partition() )},
        ghost{array::make_view<int, 1>( mesh.nodes().ghost() )},
        halo{array::make_view<int, 1>( mesh.nodes().halo() )},
        node_flags{array::make_view<int, 1>( mesh.nodes().flags() )},
        water{array::make_view<int, 1>( mesh.nodes().add(
            Field( "water", array::make_datatype<int>(), array::make_shape( mesh.nodes().size() ) ) ) )},
        //core{array::make_view<int, 1>( mesh.nodes().add(
        //    Field( "core", array::make_datatype<int>(), array::make_shape( mesh.nodes().size() ) ) ) )},
        master_glb_idx{array::make_view<gidx_t, 1>( mesh.nodes().add( Field(
            "master_global_index", array::make_datatype<gidx_t>(), array::make_shape( mesh.nodes().size() ) ) ) )} {}
};

struct Cells {
    array::ArrayView<int, 1> part;
    array::ArrayView<int, 1> halo;
    array::ArrayView<gidx_t, 1> glb_idx;
    array::ArrayView<int, 1> flags_view;
    mesh::HybridElements::Connectivity& node_connectivity;
    util::detail::BitflagsView<int> flags( idx_t i ) { return util::Topology::view( flags_view( i ) ); }
    Cells( Mesh& mesh ) :
        part{array::make_view<int, 1>( mesh.cells().partition() )},
        halo{array::make_view<int, 1>( mesh.cells().halo() )},
        glb_idx{array::make_view<gidx_t, 1>( mesh.cells().global_index() )},
        flags_view{array::make_view<int, 1>( mesh.cells().flags() )},
        node_connectivity( mesh.cells().node_connectivity() ) {}
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

void OrcaMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    ATLAS_TRACE( "OrcaMeshGenerator::generate" );
    using Topology = util::Topology;

    OrcaGrid orca{grid};
    ATLAS_ASSERT( orca );
    ATLAS_ASSERT( !mesh.generated() );


    Configuration SR_cfg;
    SR_cfg.mypart = mypart_;
    SurroundingRectangle SR( grid, distribution, SR_cfg );


    // clone some grid properties
    setGrid( mesh, grid, distribution );

    OrcaGrid rg{grid};

    int nparts     = nparts_;
    int nx         = rg.nx();
    int ny         = rg.ny();
    int ny_halo    = ny + orca.haloNorth();
    int iy_glb_max = ny + orca.haloNorth() - 1;
    int iy_glb_min = -orca.haloSouth();
    int ix_glb_max = nx + orca.haloEast() - 1;
    int ix_glb_min = -orca.haloWest();

    int nx_halo_WE = nx + orca.haloEast() + orca.haloWest();
    int ny_halo_NS = ny + orca.haloNorth() + orca.haloSouth();

    // vector of local indices: necessary for remote indices of ghost nodes
    idx_t glbarray_offset  = -( nx_halo_WE * iy_glb_min ) - ix_glb_min;
    idx_t glbarray_jstride = nx_halo_WE;

    auto index_glbarray = [&]( idx_t i, idx_t j ) {
        ATLAS_ASSERT( i <= ix_glb_max );
        ATLAS_ASSERT( j <= iy_glb_max );
        return glbarray_offset + j * glbarray_jstride + i;
    };

    std::string partitioner_name = distribution.type();
    const bool serial_distribution = (nparts == 1 || partitioner_name == "serial");

    auto partition = [&]( idx_t i, idx_t j ) -> int {
        if ( nparts == 1 ) {
            return 0;
        }
        auto clamp = []( idx_t value, idx_t lower, idx_t upper ) {
            // in C++17 this is std::clamp
            return std::max( lower, std::min( value, upper ) );
        };
        i = clamp( i, 0, nx - 1 );
        j = clamp( j, 0, ny - 1 );
        return distribution.partition( j * nx + i );
    };

    //---------------------------------------------------

    int nnodes = SR.nb_nodes;
    int ncells = SR.nb_cells;

    if ( serial_distribution ) {
        ATLAS_ASSERT( ( nx_halo_WE * ny_halo_NS ) == nnodes );
    }

    // define nodes and associated properties
    mesh.nodes().resize( nnodes );
    Nodes nodes( mesh );

    // define cells and associated properties
#if ATLAS_TEMPORARY_ELEMENTTYPES
    // DEPRECATED
    mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
#else
    // Use this since atlas 0.35.0
    mesh.cells().add( mesh::ElementType::create("Quadrilateral"), ncells );
#endif

    Cells cells( mesh );

    int inode_nonghost, inode_ghost;


    int ix_pivot = nx / 2;
    bool patch   = not orca.ghost( ix_pivot + 1, ny - 1 );

    std::vector<idx_t> node_index( SR.size, -1 );

    {
        ATLAS_TRACE( "nodes" );

        // loop over nodes and set properties
        inode_nonghost = 0;
        inode_ghost    = SR.nb_nodes_owned;  // ghost nodes start counting after nonghost nodes

        ATLAS_TRACE_SCOPE( "indexing" )
        for ( idx_t iy = 0; iy < SR.ny; iy++ ) {
            idx_t iy_glb = SR.iy_min + iy;
            ATLAS_ASSERT( iy_glb < ny_halo );
            for ( idx_t ix = 0; ix < SR.nx; ix++ ) {
                idx_t ii = SR.index( ix, iy );
                // node properties
                if ( SR.is_node[ii] ) {
                    // set node counter
                    if ( SR.is_ghost[ii] ) {
                        node_index[ii] = inode_ghost++;
                        ATLAS_ASSERT( node_index[ii] < SR.nb_nodes );
                    }
                    else {
                        node_index[ii] = inode_nonghost++;
                        ATLAS_ASSERT( node_index[ii] < SR.nb_nodes_owned );
                    }
                }
            }
        }

        inode_nonghost = 0;
        inode_ghost    = SR.nb_nodes_owned;  // ghost nodes start counting after nonghost nodes

        ATLAS_TRACE_SCOPE( "filling" )
        atlas_omp_parallel_for( idx_t iy = 0; iy < SR.ny; iy++ ) {
            idx_t iy_glb = SR.iy_min + iy;
            ATLAS_ASSERT( iy_glb < ny_halo );
            double lon00 = orca.xy( 0, 0 ).x();
            double west  = lon00 - 90.;

            auto normalise_lon00 = util::NormaliseLongitude( lon00 - 180. );
            double lon1          = normalise_lon00( orca.xy( 1, iy_glb ).x() );
            if ( lon1 < lon00 - 10. ) {
                west = lon00 - 20.;
            }

            auto normalise_lon_first_half  = util::NormaliseLongitude{west};
            auto normalise_lon_second_half = util::NormaliseLongitude{lon00 + 90.};
            for ( idx_t ix = 0; ix < SR.nx; ix++ ) {
                idx_t ix_glb   = SR.ix_min + ix;
                idx_t ix_glb_master = ix_glb;
                auto normalise = [&]( double _xy[2] ) {
                    if ( ix_glb_master < nx / 2 ) {
                        _xy[LON] = normalise_lon_first_half( _xy[LON] );
                    }
                    else {
                        _xy[LON] = normalise_lon_second_half( _xy[LON] );
                    }
                };

                idx_t ii = SR.index( ix, iy );
                // node properties
                if ( SR.is_node[ii] ) {
                    idx_t inode = node_index[ii];

                    // ghost nodes
                    nodes.ghost( inode ) = SR.is_ghost[ii];
                    if ( iy_glb > 0 or ix_glb < 0 ) {
                        nodes.ghost( inode ) = nodes.ghost( inode ) || orca.ghost( ix_glb, iy_glb );
                    }

                    // flags
                    auto flags = nodes.flags( inode );
                    flags.reset();

                    // global index
                    nodes.glb_idx( inode ) = index_glbarray( ix_glb, iy_glb ) + 1;  // no periodic point

                    // grid coordinates

                    nodes.ij( inode, XX ) = ix_glb;
                    nodes.ij( inode, YY ) = iy_glb;

                    double _xy[2];

                    // normal calculation
                    orca.xy( ix_glb, iy_glb, _xy );
                    normalise( _xy );

                    nodes.xy( inode, LON ) = _xy[LON];
                    nodes.xy( inode, LAT ) = _xy[LAT];

                    // geographic coordinates by using projection
                    nodes.lonlat( inode, LON ) = _xy[LON];
                    nodes.lonlat( inode, LAT ) = _xy[LAT];

                    // part and remote_idx
                    nodes.part( inode )           = SR.parts[ii];
                    nodes.remote_idx( inode )     = inode;
                    nodes.master_glb_idx( inode ) = nodes.glb_idx( inode );
                    if ( nodes.ghost( inode ) ) {
                        gidx_t master_idx             = orca.periodicIndex( ix_glb, iy_glb );
                        nodes.master_glb_idx( inode ) = master_idx + 1;
                        idx_t master_i, master_j;
                        orca.index2ij( master_idx, master_i, master_j );
                        nodes.part( inode ) = partition( master_i, master_j );
                        flags.set( Topology::GHOST );
                        nodes.remote_idx( inode ) = serial_distribution ? master_idx : -1;

                        if( nodes.glb_idx(inode) != nodes.master_glb_idx(inode) ) {
                            if ( ix_glb >= nx - orca.haloWest() ) {
                                flags.set( Topology::PERIODIC );
                            }
                            else if ( ix_glb < orca.haloEast() - 1 ) {
                                flags.set( Topology::PERIODIC );
                            }
                            if ( iy_glb >= ny - orca.haloNorth() - 1 ) {
                                flags.set( Topology::PERIODIC );
                                if( _xy[LON] > lon00 + 90. ) {
                                    flags.set( Topology::EAST );
                                }
                                else {
                                    flags.set( Topology::WEST );
                                }
                            }

                            if( flags.check( Topology::PERIODIC ) ) {
                                // It can still happen that nodes were flagged as periodic wrongly
                                // e.g. where the grid folds into itself

                                idx_t iy_glb_master;
                                double xy_master[2];
                                orca.index2ij( master_idx, ix_glb_master, iy_glb_master );
                                orca.lonlat(ix_glb_master,iy_glb_master,xy_master);
                                normalise( xy_master );
                                if( std::abs(xy_master[LON] - _xy[LON]) < 1.e-12 ) {
                                    flags.unset(Topology::PERIODIC);
                                }
                            }

                        }

                    }

                    flags.set( orca.land( ix_glb, iy_glb ) ? Topology::LAND : Topology::WATER );

                    if ( ix_glb <= 0 ) {
                        flags.set( Topology::BC | Topology::WEST );
                    }
                    else if ( ix_glb >= nx ) {
                        flags.set( Topology::BC | Topology::EAST );
                    }

                    nodes.water( inode ) = orca.water( ix_glb, iy_glb );
                    //nodes.core( inode ) = not orca.ghost( ix_glb, iy_glb );
                    nodes.halo( inode ) = [&]() -> int {
                        if ( ix_glb < 0 ) {
                            return -ix_glb;
                        }
                        else if ( ix_glb > nx ) {
                            return ix_glb - nx;
                        }
                        else if ( iy_glb < 0 ) {
                            return 0;
                        }
                        else if ( iy_glb >= ny ) {
                            int h = 1;
                            if ( patch && ix_glb < ix_pivot ) {  // case of eg ORCA1_T
                                h = 0;
                            }
                            return h;
                        }
                        else {
                            return 0;
                        }
                    }();
                }
            }
        }
    }
    std::vector<idx_t> cell_index( SR.size );
    // loop over nodes and define cells
    {
        ATLAS_TRACE( "elements" );
        idx_t jcell = 0;
        ATLAS_TRACE_SCOPE( "indexing" );
        for ( idx_t iy = 0; iy < SR.ny - 1; iy++ ) {      // don't loop into ghost/periodicity row
            for ( idx_t ix = 0; ix < SR.nx - 1; ix++ ) {  // don't loop into ghost/periodicity column
                idx_t ii = SR.index( ix, iy );
                if ( !SR.is_ghost[ii] ) {
                    cell_index[ii] = jcell++;
                }
            }
        }

        ATLAS_TRACE_SCOPE( "filling" )
        atlas_omp_parallel_for( idx_t iy = 0; iy < SR.ny - 1; iy++ ) {  // don't loop into ghost/periodicity row
            for ( idx_t ix = 0; ix < SR.nx - 1; ix++ ) {                // don't loop into ghost/periodicity column
                idx_t ii   = SR.index( ix, iy );
                int ix_glb = SR.ix_min + ix;
                int iy_glb = SR.iy_min + iy;
                if ( !SR.is_ghost[ii] ) {
                    idx_t jcell = cell_index[ii];

                    // define cell corners (local indices)
                    std::array<idx_t, 4> quad_nodes;
                    quad_nodes[0] = node_index[SR.index( ix, iy )];          // lower left
                    quad_nodes[1] = node_index[SR.index( ix + 1, iy )];      // lower right
                    quad_nodes[2] = node_index[SR.index( ix + 1, iy + 1 )];  // upper right
                    quad_nodes[3] = node_index[SR.index( ix, iy + 1 )];      // upper left

                    cells.flags( jcell ).reset();

                    cells.node_connectivity.set( jcell, quad_nodes.data() );
                    cells.part( jcell )    = nodes.part( quad_nodes[0] );
                    cells.glb_idx( jcell ) = ( iy_glb - iy_glb_min ) * ( nx_halo_WE - 1 ) + ( ix_glb - ix_glb_min ) + 1;

                    if ( iy_glb >= ny - 1 ) {
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
                            else if ( ix_glb >= nx ) {
                                h = ix_glb - ( nx - 1 );
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
                    if ( orca.invalidElement( SR.ix_min + ix, SR.iy_min + iy ) ) {
                        cells.flags( jcell ).set( Topology::INVALID );
                    }
                }
            }
        }
    }

    if ( serial_distribution ) {
        // Bypass for "BuildParallelFields"
        mesh.nodes().metadata().set( "parallel", true );

        // Bypass for "BuildPeriodicBoundaries"
        mesh.metadata().set( "periodic", true );
    }

    mesh.nodes().metadata().set<size_t>( "NbRealPts", nnodes );
    mesh.nodes().metadata().set<size_t>( "NbVirtualPts", size_t( 0 ) );
}

OrcaMeshGenerator::OrcaMeshGenerator( const eckit::Parametrisation& config ) {
    config.get( "partition", mypart_ = mpi::rank() );
    config.get( "partitions", nparts_ = mpi::size() );
}

void OrcaMeshGenerator::generate( const Grid& grid, const grid::Partitioner& partitioner, Mesh& mesh ) const {
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

}  // namespace meshgenerator
}  // namespace orca
}  // namespace atlas
