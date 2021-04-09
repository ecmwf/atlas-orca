
#include "OrcaMeshGenerator.h"

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
#include "atlas/interpolation/element/Quad3D.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "FixupMesh.h"


namespace atlas {
namespace orca {
namespace meshgenerator {

struct PointIJ {
    int i,j;
};

std::vector<PointIJ> detect_folds(const OrcaGrid& grid) {
    Log::info() << "Detecting folds..." << std::endl;
    std::vector<PointIJ> folds;

    geometry::UnitSphere unit_sphere;

    constexpr idx_t j0 = 0;
    double lon0 = grid->lonlat(0,j0).lon();
    auto normalise_lon_first_half  = util::NormaliseLongitude{lon0 - 90.};
    auto normalise_lon_second_half = util::NormaliseLongitude{lon0 + 90.};

    auto normalized_lon = [&](idx_t i, idx_t j) {
        double lon = grid->lonlat(i,j).lon();
        if( i < grid.nx()/2 ) {
            return normalise_lon_first_half(lon);
        }
        else {
            return normalise_lon_second_half(lon);
        }
    };//util::NormaliseLongitude(lon_0);
    using atlas::interpolation::element::Quad3D;
    for( idx_t i=grid.haloWest(); i < grid.nx()+grid.haloEast()-1; ++i ) {

        double lon_i = normalized_lon(i,j0);
        double lon_ip1 = normalized_lon(i+1,j0);
        if( lon_i > lon_ip1 ) {
            idx_t j = j0;
            idx_t j_invalid = j;
            for( j = j0; j < grid.ny(); ++j ) {
                lon_i = normalized_lon(i,j);
                lon_ip1 = normalized_lon(i+1,j);
                std::array<PointXYZ,4> p;
                p[0] = unit_sphere.xyz(grid.lonlat(i,j));
                p[1] = unit_sphere.xyz(grid.lonlat(i,j-1));
                p[2] = unit_sphere.xyz(grid.lonlat(i+1,j-1));
                p[3] = unit_sphere.xyz(grid.lonlat(i+1,j));
                if( not Quad3D{p[0],p[1],p[2],p[3]}.validate() ) {
                    j_invalid = std::max(j,j_invalid);
                }
                if( grid.water(i,j) ) {
                    break;
                }
            }
            folds.emplace_back(PointIJ{i,j_invalid});
            Log::info() << "Found fold NW={"<<folds.back().i<<","<<folds.back().j<<"}" << std::endl;
        }
    }
    return folds;
}


struct Configuration {
    bool include_south_pole;
    int nparts;
    int mypart;
};

struct SurroundingRectangle {
    std::vector<int> parts;
    std::vector<bool> is_ghost;
    std::vector<bool> is_node;
    int size;
    int nx, ny;
    int ix_min, ix_max;
    int iy_min, iy_max;
    int nb_nodes;
    int nb_cells;
    int nb_nodes_sp;
    int nb_cells_sp;
    int nb_nodes_owned;
    bool include_south_pole;

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
            std::vector<int> nb_nodes_owned_TP( atlas_omp_get_max_threads(), 0 );
            atlas_omp_parallel {
                const size_t tid = atlas_omp_get_thread_num();
                atlas_omp_for( idx_t iy = iy_glb_min; iy <= iy_glb_max; iy++ ) {
                    for ( idx_t ix = ix_glb_min; ix <= ix_glb_max; ix++ ) {
                        int p = partition( ix, iy );
                        if ( p == mypart ) {
                            ix_min = std::min( ix_min, ix );
                            ix_max = std::max( ix_max, ix );
                            iy_min = std::min( iy_min, iy );
                            iy_max = std::max( iy_max, iy );
                            nb_nodes_owned_TP[tid]++;
                        }
                    }
                }
            }
            nb_nodes_owned = std::accumulate( nb_nodes_owned_TP.begin(), nb_nodes_owned_TP.end(), 0 );
        }

        // add one row/column for ghost nodes (which include periodicity points)
        ix_max = std::min( ix_glb_max, ix_max + 1 );
        iy_max = std::min( iy_glb_max, iy_max + 1 );

        // dimensions of surrounding rectangle (SR)
        nx = ix_max - ix_min + 1;
        ny = iy_max - iy_min + 1;

        // upper estimate for number of nodes
        include_south_pole = cfg.include_south_pole && iy_min == iy_glb_min;
        size               = nx * ( ny + ( include_south_pole ? 1 : 0 ) );

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
            if ( include_south_pole ) {
                int ii_south = nx * ny;
                for ( idx_t ix = 0; ix < nx; ix++ ) {
                    is_ghost.at( ii_south + ix ) = false;
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
            nb_cells    = 0;
            nb_nodes    = 0;
            nb_cells_sp = 0;
            nb_nodes_sp = 0;
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
            if ( include_south_pole ) {
                for ( idx_t ix = 0; ix < nx; ix++ ) {
                    idx_t ii       = index_sp( ix );
                    idx_t ii_above = index( ix, 0 );
                    is_node[ii]    = is_node[ii_above];
                    is_ghost[ii]   = is_ghost[ii_above];
                    if ( is_node[ii] ) {
                        ++nb_nodes_sp;
                    }
                }
                nb_cells_sp = ( nb_nodes_sp > 1 ) ? nb_nodes_sp - 1 : 0;
                nb_cells += nb_cells_sp;
                nb_nodes += nb_nodes_sp;
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
            "master_global_index", array::make_datatype<gidx_t>(), array::make_shape( mesh.nodes().size() ) ) ) )}
    {}
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
    ATLAS_TRACE( "OrcaMeshGenerator::generate" );
    using Topology = util::Topology;

    OrcaGrid orca{grid};
    ATLAS_ASSERT( orca );
    ATLAS_ASSERT( !mesh.generated() );

    auto folds = detect_folds(orca);

    Configuration SR_cfg;
    SR_cfg.include_south_pole = include_pole_;
    SR_cfg.mypart             = mypart_;
    SR_cfg.nparts             = nparts_;
    SurroundingRectangle SR( grid, distribution, SR_cfg );


    // clone some grid properties
    setGrid( mesh, grid, distribution );

    OrcaGrid rg{grid};

    bool include_south_pole = include_pole_;

    int mypart     = mypart_;
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

    int nnodes         = SR.nb_nodes;
    int ncells         = SR.nb_cells;
    int nb_extra_nodes = SR.nb_nodes_sp;  // By South Pole
    int nb_extra_cells = SR.nb_cells_sp;  // By South Pole

    if ( nparts == 1 ) {
        ATLAS_ASSERT( ( nx_halo_WE * ny_halo_NS ) + nb_extra_nodes == nnodes );
    }

    // define nodes and associated properties
    mesh.nodes().resize( nnodes );
    Nodes nodes( mesh );

    // define cells and associated properties
    mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
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
            double lon0                    = orca.xy( 0, iy_glb ).x();
            auto normalise_lon_first_half  = util::NormaliseLongitude{lon0 - 90.};
            auto normalise_lon_second_half = util::NormaliseLongitude{lon0 + 90.};
            for ( idx_t ix = 0; ix < SR.nx; ix++ ) {
                idx_t ix_glb   = SR.ix_min + ix;
                auto normalise = [&]( double _xy[2] ) {
                    if ( ix_glb < nx / 2 ) {
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
                        nodes.remote_idx( inode ) = -1;
                    }

                    if ( ix_glb >= nx - orca.haloWest() ) {
                        flags.set( Topology::PERIODIC );
                    }
                    else if ( ix_glb < orca.haloEast() ) {
                        flags.set( Topology::PERIODIC );
                    }
                    if ( iy_glb >= ny - orca.haloNorth() - 1 ) {
                        flags.set( Topology::PERIODIC );
                    }

                    flags.set( orca.land( ix_glb, iy_glb ) ? Topology::LAND : Topology::WATER );

                    if ( ix_glb <= 0 ) {
                        flags.set( Topology::BC | Topology::WEST );
                    }
                    else if ( ix_glb >= nx ) {
                        flags.set( Topology::BC | Topology::EAST );
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
            }
        }
    }
    idx_t inode_south = nnodes - nb_extra_nodes;
    if ( SR.include_south_pole ) {
        idx_t inode                    = inode_south;
        gidx_t glb_idx_0               = index_glbarray( ix_glb_max, iy_glb_max ) + 2;
        double lon0                    = orca.xy( 0, -orca.haloSouth() ).x();
        auto normalise_lon_first_half  = util::NormaliseLongitude{lon0 - 90.};
        auto normalise_lon_second_half = util::NormaliseLongitude{lon0 + 90.};

        for ( idx_t ix = 0; ix < SR.nx; ix++ ) {
            idx_t ix_glb   = SR.ix_min + ix;
            auto normalise = [&]( double _xy[2] ) {
                if ( ix_glb < nx / 2 ) {
                    _xy[LON] = normalise_lon_first_half( _xy[LON] );
                }
                else {
                    _xy[LON] = normalise_lon_second_half( _xy[LON] );
                }
            };


            idx_t ii = SR.index_sp( ix );
            if ( ii >= SR.is_node.size() ) {
                throw_OutOfRange( "ii", ii, SR.is_node.size(), Here() );
            }

            if ( SR.is_node[ii] ) {
                node_index[ii]    = inode;
                idx_t inode_north = node_index[SR.index( ix, 0 )];

                PointXY p{orca.xy( ix_glb, -orca.haloSouth() ).x(), -90.};
                normalise( p.data() );

                if ( inode >= nnodes ) {
                    throw_OutOfRange( "inode", inode, nnodes, Here() );
                }

                nodes.glb_idx( inode ) = glb_idx_0 + index_glbarray( ix_glb, iy_glb_min );

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
        }
    }

    idx_t icell_south = 0;
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
        icell_south = jcell;

        ATLAS_TRACE_SCOPE( "filling" )
        atlas_omp_parallel_for( idx_t iy = 0; iy < SR.ny - 1; iy++ ) {  // don't loop into ghost/periodicity row
            for ( idx_t ix = 0; ix < SR.nx - 1; ix++ ) {                // don't loop into ghost/periodicity column
                idx_t ii   = SR.index( ix, iy );
                int ix_glb = SR.ix_min + ix;
                int iy_glb = SR.iy_min + iy;
                if ( !SR.is_ghost[ii] ) {
                    idx_t jcell = cell_index[ii];

                    // define cell corners (local indices)
                    std::array<int, 4> quad_nodes;
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
                }
            }
        }
    }

    // Special modifications
    {
        ATLAS_TRACE( "pole elements" );
        if ( SR.include_south_pole ) {
            gidx_t glb_idx_max{0};
            idx_t inode = inode_south;
            idx_t jcell = icell_south;
            for ( idx_t ix = 0; ix < SR.nx - 1; ++ix ) {
                idx_t ix_glb = SR.ix_min + ix;
                if ( !SR.is_ghost[SR.index( ix, 0 )] ) {
                    ATLAS_ASSERT( jcell < mesh.cells().size() );
                    std::array<int, 4> quad_nodes;
                    quad_nodes[0] = node_index[SR.index( ix, 0 )];
                    quad_nodes[3] = node_index[SR.index( ix + 1, 0 )];
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
        if( not folds.empty() ) {
            for( auto& ij: folds ) {
                idx_t ix_fold = ij.i - SR.ix_min;
                idx_t iy_fold = ij.j - SR.iy_min;
                if( ix_fold >= 0 && ix_fold < SR.nx && iy_fold >= 0 && iy_fold < SR.ny ) {
                    for( idx_t iy = 0; iy <= iy_fold; ++iy ) {
                        idx_t ii   = SR.index( ix_fold, iy );
                        if ( !SR.is_ghost[ii] ) {
                            idx_t jcell = cell_index[ii];
                            cells.flags( jcell ).set( Topology::INVALID );
                        }
                    }
                }
            }
        }
        if ( fixup_ ) {
            FixupMesh::create( util::Config( "type", grid.name() ) )->execute( mesh );
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
