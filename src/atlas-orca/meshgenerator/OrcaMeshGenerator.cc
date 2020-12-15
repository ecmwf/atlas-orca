
#include "OrcaMeshGenerator.h"

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
#include "atlas/meshgenerator/detail/RegularMeshGenerator.h"
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

#include "FixupMesh.h"
#include "FixupMeshUtils.h"

#define DEBUG_OUTPUT_DETAIL 0
#define DEBUG_OUTPUT 0

namespace atlas {
namespace orca {
namespace meshgenerator {

struct SurroundingRectangle {
    std::vector<int> local_idx;
    std::vector<int> parts;
    std::vector<bool> is_ghost;
    std::vector<bool> is_node;
    int nnodes;
    int nx, ny;
    int ix_min, ix_max;
    int iy_min, iy_max;
    bool include_south_pole;
    int index( int i, int j ) { return j * nx + i; }
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
    array::ArrayView<int, 1> lsm;
    array::ArrayView<int, 1> core;
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
        lsm{array::make_view<int, 1>( mesh.nodes().add(
            Field( "lsm", array::make_datatype<int>(), array::make_shape( mesh.nodes().size() ) ) ) )},
        core{array::make_view<int, 1>( mesh.nodes().add(
            Field( "core", array::make_datatype<int>(), array::make_shape( mesh.nodes().size() ) ) ) )},
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

    // TODO: This hard-coding should instead go to the data-file, and become queried via the OrcaGrid
    using t = std::tuple<std::string, bool>;
    std::vector<t> resolutions{t{"025", false}, t{"12", false}, t{"1", true}, t{"2", false}};

    std::vector<std::string> prefixes = {"ORCA", "eORCA"};
    std::map<std::string, bool> patch;
    for ( auto& resol : resolutions ) {
        std::string resolution;
        bool patch_T;
        std::tie( resolution, patch_T ) = resol;
        for ( auto& prefix : prefixes ) {
            patch[prefix + resolution + "_U"] = patch_T;
            patch[prefix + resolution + "_V"] = not patch_T;
            patch[prefix + resolution + "_W"] = patch_T;
            patch[prefix + resolution + "_T"] = patch_T;
            patch[prefix + resolution + "_F"] = not patch_T;
        }
    }

    // Mimic hole in South pole, and numbering from South to North. patch determines if endpoint is at North Pole
    StructuredGrid::YSpace yspace{grid::LinearSpacing{{-80., 90.}, orca.ny(), not patch.at( orca.name() )}};
    // Periodic xspace
    StructuredGrid::XSpace xspace{grid::LinearSpacing{{0., 360.}, orca.nx(), false}};

    return StructuredGrid{xspace, yspace};
}
}  // namespace

void OrcaMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    ATLAS_TRACE();
    using Topology = util::Topology;

    OrcaGrid orca{grid};
    ATLAS_ASSERT( orca );
    ATLAS_ASSERT( !mesh.generated() );
    // clone some grid properties
    setGrid( mesh, grid, distribution );

    OrcaGrid rg{grid};
    SurroundingRectangle SR;

    bool include_south_pole = include_pole_;

    int mypart     = mypart_;
    int nparts     = nparts_;
    int nx         = rg.nx();
    int ny         = rg.ny();
    int ny_halo    = ny + orca.haloNorth();
    int ny_halo_NS = ny + orca.haloNorth() + orca.haloSouth();
    int iy_glb_max = ny + orca.haloNorth() - 1;
    int iy_glb_min = -orca.haloSouth();
    int ix_glb_max = nx + orca.haloEast() - 1;
    int ix_glb_min = -orca.haloWest();
    int nx_halo_WE = nx + orca.haloEast() + orca.haloWest();

    if ( nparts > 1 ) {
        //include_south_pole = false;
    }

    int nb_extra_nodes = 0;
    int ncells;

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<idx_t> local_idx( nx_halo_WE * ( ny_halo_NS + ( include_south_pole ? 1 : 0 ) ), -1 );
    idx_t local_idx_offset  = -( nx_halo_WE * iy_glb_min ) - ix_glb_min;
    idx_t local_idx_jstride = nx_halo_WE;

    auto compute_global_offset = [&]( idx_t i, idx_t j ) {
        ATLAS_ASSERT( i <= ix_glb_max );
        ATLAS_ASSERT( j <= iy_glb_max );
        return local_idx_offset + j * local_idx_jstride + i;
    };

    auto global_to_local = [&]( idx_t i, idx_t j ) {
        idx_t n = compute_global_offset( i, j );
        return local_idx.at( n );
    };


    std::vector<int> current_idx( nparts, 0 );  // index counter for each proc

    auto clamp = []( idx_t value, idx_t lower, idx_t upper ) {
        // in C++17 this is std::clamp
        return std::max( lower, std::min( value, upper ) );
    };
    auto partition = [&]( idx_t i, idx_t j ) -> int {
        i = clamp( i, 0, nx - 1 );
        j = clamp( j, 0, ny - 1 );
        return distribution.partition( j * nx + i );
    };

    // determine rectangle (ix_min:ix_max) x (iy_min:iy_max) surrounding the nodes on this processor
    idx_t ix_glb, iy_glb, ix, iy;
    int nnodes_nonghost, nnodes;  // number of nodes: non-ghost; total;  inside surrounding rectangle
    int ii;
    int ii_glb_south;
    // loop over all points to determine local indices and surroundig rectangle
    SR.ix_min       = nx + 1;
    SR.ix_max       = 0;
    SR.iy_min       = ny + 1;
    SR.iy_max       = 0;
    nnodes_nonghost = 0;

    {
        gidx_t ii_glb = 0;
        for ( iy = iy_glb_min; iy <= iy_glb_max; iy++ ) {
            for ( ix = ix_glb_min; ix <= ix_glb_max; ix++ ) {
                int p                  = partition( ix, iy );
                local_idx.at( ii_glb ) = current_idx[p]++;  // store local index on the local proc of this point
                if ( p == mypart ) {
                    ++nnodes_nonghost;  // non-ghost node: belongs to this part
                    SR.ix_min = std::min( SR.ix_min, ix );
                    SR.ix_max = std::max( SR.ix_max, ix );
                    SR.iy_min = std::min( SR.iy_min, iy );
                    SR.iy_max = std::max( SR.iy_max, iy );
                }
                ++ii_glb;  // global index
            }
        }
        ii_glb_south = ii_glb;
        if ( include_south_pole ) {
            iy = iy_glb_min - 1;
            for ( ix = ix_glb_min; ix <= ix_glb_max; ix++ ) {
                int p                  = partition( ix, iy + 1 );
                local_idx.at( ii_glb ) = current_idx[p]++;  // store local index on the local proc of this point
                ++ii_glb;                                   // global index
            }
        }
    }

    // add one row/column for ghost nodes (which include periodicity points)
    SR.ix_max = std::min( ix_glb_max, SR.ix_max + 1 );
    SR.iy_max = std::min( iy_glb_max, SR.iy_max + 1 );

    // dimensions of surrounding rectangle (SR)
    //int nxl = ix_max - ix_min + 1;
    //int nyl = iy_max - iy_min + 1;
    SR.nx = SR.ix_max - SR.ix_min + 1;
    SR.ny = SR.iy_max - SR.iy_min + 1;

    // upper estimate for number of nodes
    SR.include_south_pole = include_south_pole && SR.iy_min == iy_glb_min;
    SR.nnodes             = SR.nx * ( SR.ny + ( SR.include_south_pole ? 1 : 0 ) );

    // partitions and local indices in SR
    SR.parts.resize( SR.nnodes, -1 );
    SR.local_idx.resize( SR.nnodes, -1 );
    SR.is_ghost.resize( SR.nnodes, true );
    // vectors marking nodes that are necessary for this proc's cells
    SR.is_node.resize( SR.nnodes, false );

    ii = 0;  // index inside SR
    for ( iy = 0; iy < SR.ny; iy++ ) {
        iy_glb = SR.iy_min + iy;  // global y-index
        for ( ix = 0; ix < SR.nx; ix++ ) {
            ix_glb            = SR.ix_min + ix;  // global x-index
            SR.parts.at( ii ) = partition( ix_glb, iy_glb );

            if ( iy_glb < ny_halo ) {
                SR.local_idx.at( ii ) = global_to_local( ix_glb, iy_glb );
                SR.is_ghost.at( ii )  = ( SR.parts.at( ii ) != mypart );
            }
            ++ii;
        }
    }
    int ii_south = ii;
    if ( SR.include_south_pole ) {
        for ( ix = 0; ix < SR.nx; ix++ ) {
            ix_glb                = SR.ix_min + ix;
            SR.local_idx.at( ii ) = local_idx.at( ii_glb_south + ix_glb );
            SR.is_ghost.at( ii )  = false;
            ++ii;
        }
    }
    ATLAS_ASSERT( ii == SR.nx * ( SR.ny + ( SR.include_south_pole ? 1 : 0 ) ) );

    // determine number of cells and number of nodes
    nnodes = 0;
    ncells = 0;
    for ( iy = 0; iy < SR.ny - 1; iy++ ) {  // don't loop into ghost/periodicity row
        iy_glb = SR.iy_min + iy;
        for ( ix = 0; ix < SR.nx - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ix_glb = SR.ix_min + ix;
            ii     = SR.index( ix, iy );
            if ( !SR.is_ghost[ii] ) {
                ++ncells;

                // mark this node as being used
                if ( !SR.is_node[ii] ) {
                    ++nnodes;
                    SR.is_node[ii] = true;
                }
                // mark lowerright corner
                ii = SR.index( ix + 1, iy );
                if ( !SR.is_node[ii] ) {
                    ++nnodes;
                    SR.is_node[ii] = true;
                }
                // mark upperleft corner
                ii = SR.index( ix, iy + 1 );
                if ( !SR.is_node[ii] ) {
                    ++nnodes;
                    SR.is_node[ii] = true;
                }
                // mark upperright corner
                ii = SR.index( ix + 1, iy + 1 );
                if ( !SR.is_node[ii] ) {
                    ++nnodes;
                    SR.is_node[ii] = true;
                }
            }
        }
    }
    nb_extra_nodes = 0;
    if ( SR.include_south_pole ) {
        ii = ii_south;
        for ( ix = 0; ix < SR.nx; ix++ ) {
            ix_glb             = SR.ix_min + ix;
            idx_t ii_row_north = SR.index( ix, 0 );
            SR.is_node[ii]     = SR.is_node[ii_row_north];
            SR.is_ghost[ii]    = SR.is_ghost[ii_row_north];
            if ( SR.is_node[ii] ) {
                ++nb_extra_nodes;
            }
            ++ii;
        }
    }
    int nb_extra_cells = ( nb_extra_nodes > 1 ) ? nb_extra_nodes - 1 : 0;

    nnodes += nb_extra_nodes;
    ncells += nb_extra_cells;
    if ( nparts == 1 ) {
        ATLAS_ASSERT( ( nx_halo_WE * ny_halo_NS ) + nb_extra_nodes == nnodes );
    }

    // define nodes and associated properties
    mesh.nodes().resize( nnodes );
    Nodes nodes( mesh );

    // define cells and associated properties
    mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
    int quad_begin = mesh.cells().elements( 0 ).begin();
    Cells cells( mesh );

    std::array<int, 4> quad_nodes;
    int jcell = quad_begin;
    int inode, inode_nonghost, inode_ghost;

    // loop over nodes and set properties
    ii             = 0;
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes

    int ix_pivot = nx / 2;
    bool patch   = not orca.ghost( ix_pivot + 1, ny - 1 );

    struct PointIJ {
        idx_t i;
        idx_t j;
    };
    auto fold = [&]( idx_t i, idx_t j ) {
        PointIJ p;
        if ( not patch && j >= ny - 1 ) {  // ORCA2_T, ORCA025_T, ORCA12_T
            int iy_pivot = ny - 1;
            p.i          = ix_pivot - ( i - ix_pivot );
            p.j          = iy_pivot - ( j - iy_pivot );
        }
        else if ( patch && j >= ny ) {  // ORCA1_T
            double iy_pivot = double( ny ) - 0.5;
            p.i             = 2 * ix_pivot - 1 - i;
            p.j             = int( 2. * iy_pivot ) - j;
        }
        else {
            p.i = i;
            p.j = j;
        }
        return p;
    };

    for ( iy = 0; iy < SR.ny; iy++ ) {
        iy_glb = SR.iy_min + iy;
        ATLAS_ASSERT( iy_glb < ny_halo );
        double lon0                    = orca.xy( 0, iy_glb ).x();
        auto normalise_lon_first_half  = util::NormaliseLongitude{lon0 - 90.};
        auto normalise_lon_second_half = util::NormaliseLongitude{lon0 + 90.};
        auto normalise                 = [&]( double _xy[2] ) {
            if ( ix_glb < nx / 2 ) {
                _xy[LON] = normalise_lon_first_half( _xy[LON] );
            }
            else {
                _xy[LON] = normalise_lon_second_half( _xy[LON] );
            }
        };
        for ( ix = 0; ix < SR.nx; ix++ ) {
            ix_glb = ( SR.ix_min + ix );  // don't take modulus here: periodicity points
                                          // have their own global index.
            // node properties
            if ( SR.is_node[ii] ) {
                // set node counter
                if ( SR.is_ghost[ii] ) {
                    inode = inode_ghost++;
                    ATLAS_ASSERT( inode < nnodes );
                }
                else {
                    inode = inode_nonghost++;
                    ATLAS_ASSERT( inode < nnodes_nonghost );
                }

                // ghost nodes
                nodes.ghost( inode ) = SR.is_ghost[ii];
                if ( iy_glb > 0 or ix_glb < 0 ) {
                    nodes.ghost( inode ) = nodes.ghost( inode ) || orca.ghost( ix_glb, iy_glb );
                }

                // flags
                nodes.flags( inode ).reset();

                // global index
                nodes.glb_idx( inode ) = compute_global_offset( ix_glb, iy_glb ) + 1;  // no periodic point

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
                    gidx_t master_idx             = orca->periodicIndex( ix_glb, iy_glb );
                    nodes.master_glb_idx( inode ) = master_idx + 1;
                    PointIJ master;
                    orca->index2ij( master_idx, master.i, master.j );
                    nodes.part( inode ) = partition( master.i, master.j );
                    //                    ATLAS_DEBUG( "ghost: " << ix_glb << ", " << iy_glb << "  inode = " << inode
                    //                                           << " P = " << PointLonLat( _xy ) );
                    nodes.flags( inode ).set( Topology::GHOST );
                    nodes.remote_idx( inode ) = SR.local_idx[ii];
                    // change local index -- required for cells
                    SR.local_idx[ii] = inode;
                }

                if ( ix_glb >= nx - orca.haloWest() ) {
                    nodes.flags( inode ).set( Topology::PERIODIC );
                }
                else if ( ix_glb < orca.haloEast() ) {
                    nodes.flags( inode ).set( Topology::PERIODIC );
                }
                if ( iy_glb >= ny - orca.haloNorth() - 1 ) {
                    nodes.flags( inode ).set( Topology::PERIODIC );
                }

                if ( orca.land( ix_glb, iy_glb ) == 0 ) {
                    nodes.flags( inode ).set( Topology::LAND );
                }
                else {
                    nodes.flags( inode ).set( Topology::WATER );
                }
                if ( ix_glb <= 0 ) {
                    nodes.flags( inode ).set( Topology::BC | Topology::WEST );
                }
                else if ( ix_glb >= nx ) {
                    nodes.flags( inode ).set( Topology::BC | Topology::EAST );
                }

                nodes.lsm( inode )  = orca.water( ix_glb, iy_glb );
                nodes.core( inode ) = not orca.ghost( ix_glb, iy_glb );
                nodes.halo( inode ) = [&]() -> int {
                    if ( ix_glb < 0 ) {
                        return -ix_glb;
                    }
                    else if ( ix_glb > nx + 1 ) {
                        return ix_glb - nx + 1;
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
            ++ii;
        }
    }
    int inode_south = nnodes - nb_extra_nodes;
    inode           = inode_south;
    ii              = ii_south;
    if ( SR.include_south_pole ) {
        gidx_t glb_idx_0               = compute_global_offset( ix_glb_max, iy_glb_max ) + 2;
        double lon0                    = orca.xy( 0, -orca.haloSouth() ).x();
        auto normalise_lon_first_half  = util::NormaliseLongitude{lon0 - 90.};
        auto normalise_lon_second_half = util::NormaliseLongitude{lon0 + 90.};
        auto normalise                 = [&]( double _xy[2] ) {
            if ( ix_glb < nx / 2 ) {
                _xy[LON] = normalise_lon_first_half( _xy[LON] );
            }
            else {
                _xy[LON] = normalise_lon_second_half( _xy[LON] );
            }
        };

        for ( ix = 0; ix < SR.nx; ix++ ) {
            if ( ii >= SR.is_node.size() ) {
                throw_OutOfRange( "ii", ii, SR.is_node.size(), Here() );
            }

            if ( SR.is_node[ii] ) {
                idx_t inode_north = SR.local_idx[SR.index( ix, 0 )];

                ix_glb = SR.ix_min + ix;
                PointXY p{orca.xy( ix_glb, -orca.haloSouth() ).x(), -90.};
                normalise( p.data() );

                if ( inode >= nnodes ) {
                    throw_OutOfRange( "inode", inode, nnodes, Here() );
                }

                nodes.glb_idx( inode ) = glb_idx_0 + compute_global_offset( ix_glb, iy_glb_min );

                nodes.master_glb_idx( inode ) = nodes.glb_idx( inode );
                nodes.part( inode )           = nodes.part( inode_north );

                nodes.lonlat( inode, LON ) = p.x();
                nodes.lonlat( inode, LAT ) = p.y();
                nodes.xy( inode, XX )      = p.x();
                nodes.xy( inode, YY )      = p.y();
                nodes.ij( inode, XX )      = nodes.ij( inode_north, XX );
                nodes.ij( inode, YY )      = nodes.ij( inode_north, YY ) - 1;
                nodes.ghost( inode )       = nodes.ghost( inode_north );
                nodes.lsm( inode )         = nodes.lsm( inode_north );  // land
                nodes.node_flags( inode )  = nodes.node_flags( inode_north );
                nodes.halo( inode )        = nodes.halo( inode_north );
                nodes.core( inode )        = nodes.core( inode_north );
                nodes.remote_idx( inode )  = inode;
                if ( ix_glb < 0 ) {
                    nodes.master_glb_idx( inode ) = glb_idx_0 + ix_glb + nx;
                    nodes.remote_idx( inode )     = -1;
                }
                else if ( ix_glb >= nx ) {
                    nodes.master_glb_idx( inode ) = glb_idx_0 + ix_glb - nx;
                    nodes.remote_idx( inode )     = -1;
                }
                ++inode;
            }

            ++ii;
        }
    }


    // loop over nodes and define cells
    for ( iy = 0; iy < SR.ny - 1; iy++ ) {      // don't loop into ghost/periodicity row
        for ( ix = 0; ix < SR.nx - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ii         = SR.index( ix, iy );
            int ix_glb = SR.ix_min + ix;
            int iy_glb = SR.iy_min + iy;
            if ( !SR.is_ghost[ii] ) {
                // define cell corners (local indices)
                quad_nodes[0] = SR.local_idx[SR.index( ix, iy )];          // lower left
                quad_nodes[1] = SR.local_idx[SR.index( ix + 1, iy )];      // lower right
                quad_nodes[2] = SR.local_idx[SR.index( ix + 1, iy + 1 )];  // upper right
                quad_nodes[3] = SR.local_idx[SR.index( ix, iy + 1 )];      // upper left

                cells.flags( jcell ).reset();
                if ( iy + 1 >= ny ) {
                    // std::swap( quad_nodes[2], quad_nodes[3] );
                }

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

                ++jcell;
            }
        }
    }

    // Special modifications
    {
        MeshBuilder mesh_builder( mesh );
        if ( SR.include_south_pole ) {
            gidx_t glb_idx_max{0};
            inode = inode_south;
            for ( int ix = 0; ix < SR.nx - 1; ++ix ) {
                ix_glb = SR.ix_min + ix;
                if ( !SR.is_ghost[SR.index( ix, 0 )] ) {
                    ATLAS_ASSERT( jcell < mesh.cells().size() );
                    quad_nodes[0] = SR.local_idx[SR.index( ix, 0 )];
                    quad_nodes[3] = SR.local_idx[SR.index( ix + 1, 0 )];
                    quad_nodes[1] = inode;
                    quad_nodes[2] = inode + 1;
                    cells.node_connectivity.set( jcell, quad_nodes.data() );
                    cells.flags( jcell ).reset();
                    cells.flags( jcell ).set( Topology::LAND | Topology::BC | Topology::SOUTH );
                    cells.halo( jcell )    = 0;
                    cells.part( jcell )    = nodes.part( quad_nodes[0] );
                    cells.glb_idx( jcell ) = ( ny_halo_NS - 1 ) * ( nx_halo_WE - 1 ) + ( ix_glb - ix_glb_min ) + 1;
                    if ( ix_glb < 0 || ix_glb >= nx ) {
                        cells.flags( jcell ).set( Topology::GHOST );
                        cells.halo( jcell ) = ( ix_glb < 0 ? -ix_glb : ix_glb - ( nx - 1 ) );
                    }
                    glb_idx_max = std::max( glb_idx_max, cells.glb_idx( jcell ) );
                    ++jcell;
                    ++inode;
                }
            }
            ATLAS_ASSERT( jcell == mesh.cells().size() );
            if ( nparts == 1 ) {
                ATLAS_ASSERT( jcell == glb_idx_max );
            }
            mesh.metadata().set( "includes_south_pole", true );
        }
        if ( fixup_ ) {
            FixupMesh::create( util::Config( "type", grid.name() ) )->execute( mesh_builder );
        }
    }
    //nb_extra_nodes = 0;
    mesh.nodes().metadata().set( "parallel", true );
    mesh.nodes().metadata().set<size_t>( "NbRealPts", nnodes - nb_extra_nodes );
    mesh.nodes().metadata().set<size_t>( "NbVirtualPts", size_t( nb_extra_nodes ) );
}

OrcaMeshGenerator::OrcaMeshGenerator( const eckit::Parametrisation& config ) {
    bool patch_pole = false;
    config.get( "patch_pole", patch_pole );
    config.get( "include_pole", include_pole_ );
    config.get( "force_include_south_pole", include_pole_ );
    config.get( "part", mypart_ = mpi::rank() );
    config.get( "nb_parts", nparts_ = mpi::size() );
    config.get( "fixup", fixup_ );

    // This is a temporary hack as patch_pole is not implemented
    include_pole_ = include_pole_ || patch_pole;
    if ( nparts_ > 1 ) {
        fixup_ = false;  // not yet implemented
    }
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
