
#include "OrcaMeshGenerator.h"
#include "OrcaGrid.h"

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

#define DEBUG_OUTPUT_DETAIL 0
#define DEBUG_OUTPUT 0

namespace atlas {
namespace meshgenerator {

namespace {
StructuredGrid equivalent_regular_grid( const OrcaGrid& orca ) {
    ATLAS_ASSERT( orca );

    // TODO: This hard-coding should instead go to the data-file, and become queried via the OrcaGrid
    std::map<std::string, bool> patch;
    std::vector<std::tuple<std::string, bool>> resolutions = {{"025", false}, {"12", false}, {"1", true}, {"2", false}};
    std::vector<std::string> prefixes                      = {"ORCA", "eORCA"};
    for ( auto& resol : resolutions ) {
        std::string resolution;
        bool patch_T;
        std::tie( resolution, patch_T ) = resol;
        for ( auto& prefix : prefixes ) {
            patch[prefix + resolution + "_T"] = patch_T;
            patch[prefix + resolution + "_U"] = patch_T;
            patch[prefix + resolution + "_V"] = not patch_T;
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
    using Topology = util::Topology;

    OrcaGrid orca{grid};
    ATLAS_ASSERT( orca );
    ATLAS_ASSERT( !mesh.generated() );
    // clone some grid properties
    setGrid( mesh, grid, distribution );

    OrcaGrid rg{grid};

//        bool periodic_x = true;  // options.get<bool>( "periodic_x" ) or rg.periodic();
    bool periodic_x = false;  // options.get<bool>( "periodic_x" ) or rg.periodic();

    bool include_south_pole = include_pole_ ;

    int mypart  = mypart_;
    int nparts  = nparts_;
    int nx      = rg.nx();
    int ny      = rg.ny();
    int ny_halo = ny + orca->haloNorth();
    int ny_halo_NS = ny + orca->haloNorth() + orca->haloSouth();
    int iy_glb_max = ny+orca->haloNorth()-1;
    int iy_glb_min = -orca->haloSouth();
    int ix_glb_max = periodic_x ? nx-1 : nx+orca->haloEast()-1;
    int ix_glb_min = periodic_x ? 0  : -orca->haloWest();
    int nx_halo_WE = nx+orca->haloEast() + orca->haloWest();

    if( nparts > 1 ) {
        include_south_pole = false;
    }
    int nb_extra_nodes = include_south_pole ? 1 : 0;

    gidx_t ii_glb;  // global index
    int ncells;

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<idx_t> local_idx( nx_halo_WE * ny_halo_NS, -1 );
    idx_t local_idx_offset = periodic_x ? - ( nx* iy_glb_min ) : - (nx_halo_WE*iy_glb_min) - ix_glb_min;
    idx_t local_idx_jstride = (periodic_x?nx:nx_halo_WE);
    auto global_to_local = [&](idx_t i, idx_t j) {
        ATLAS_ASSERT(i <= ix_glb_max );
        ATLAS_ASSERT(j <= iy_glb_max );
        idx_t n = local_idx_offset + j*local_idx_jstride + i;
        return local_idx.at(n);
    };

    std::vector<int> current_idx( nparts, 0 );  // index counter for each proc

    auto clamp = [](idx_t value, idx_t lower, idx_t upper) {
        // in C++17 this is std::clamp
        return std::max(lower,std::min(value,upper));
    };
    auto partition = [&]( idx_t i, idx_t j ) -> int {
        i = clamp(i,0,nx-1);
        j = clamp(j,0,ny-1);
        return distribution.partition( j * nx + i );
    };

    // determine rectangle (ix_min:ix_max) x (iy_min:iy_max) surrounding the nodes
    // on this processor
    idx_t ix_min, ix_max, iy_min, iy_max, ix_glb, iy_glb, ix, iy;
    int nnodes_nonghost, nnodes;  // number of nodes: non-ghost; total;  inside
                                  // surrounding rectangle
    int nnodes_SR, ii;

    // loop over all points to determine local indices and surroundig rectangle
    ix_min          = nx + 1;
    ix_max          = 0;
    iy_min          = ny + 1;
    iy_max          = 0;
    nnodes_nonghost = 0;

    ii_glb = 0;
    for ( iy = iy_glb_min; iy <= iy_glb_max; iy++ ) {
        for ( ix = ix_glb_min; ix <= ix_glb_max; ix++ ) {
            int p             = partition( ix, iy );
            local_idx.at(ii_glb) = current_idx[p]++;  // store local index on
                                                   // the local proc of
                                                   // this point
            if ( p == mypart ) {
                ++nnodes_nonghost;  // non-ghost node: belongs to this part
                ix_min = std::min( ix_min, ix );
                ix_max = std::max( ix_max, ix );
                iy_min = std::min( iy_min, iy );
                iy_max = std::max( iy_max, iy );
            }
            ++ii_glb;  // global index
        }
    }

    // add one row/column for ghost nodes (which include periodicity points)
    if( periodic_x ) {
        ix_max++;
    }
    else {
        ix_max = std::min( ix_glb_max, ix_max + 1 );
    }
    iy_max = std::min( iy_glb_max, iy_max + 1 );

    // dimensions of surrounding rectangle (SR)
    int nxl = ix_max - ix_min + 1;
    int nyl = iy_max - iy_min + 1;

    // upper estimate for number of nodes
    nnodes_SR = nxl * nyl + nb_extra_nodes;

    // partitions and local indices in SR
    std::vector<int> parts_SR( nnodes_SR, -1 );
    std::vector<int> local_idx_SR( nnodes_SR, -1 );
    std::vector<int> is_ghost_SR( nnodes_SR, 1 );
    ii = 0;  // index inside SR
    for ( iy = 0; iy < nyl; iy++ ) {
        iy_glb = iy_min + iy;  // global y-index
        for ( ix = 0; ix < nxl; ix++ ) {
            ix_glb       = ix_min + ix;  // global x-index
            parts_SR.at(ii) = partition( ix_glb, iy_glb );

            if ( (ix_glb < nx || not periodic_x) && iy_glb < ny_halo ) {
                local_idx_SR.at(ii) = global_to_local( ix_glb, iy_glb );
                is_ghost_SR.at(ii)  = ( parts_SR.at(ii) != mypart );
            }
            else {
                local_idx_SR.at(ii) = -1;
                is_ghost_SR.at(ii)  = true;
                //                ATLAS_DEBUG( "ghost: " << ix_glb << "," << iy_glb );
            }
            ++ii;
        }
    }
    if( include_south_pole ) {
        local_idx_SR.at(ii) = -1;
        is_ghost_SR.at(ii) = true;
    }

    // vectors marking nodes that are necessary for this proc's cells
    std::vector<bool> is_node_SR( nnodes_SR, false );

    // determine number of cells and number of nodes
    nnodes     = 0;
    ncells     = 0;
    auto index = [&]( idx_t i, idx_t j ) { return j * nxl + i; };
    int nnodes_prev = 0;
    for ( iy = 0; iy < nyl - 1; iy++ ) {  // don't loop into ghost/periodicity row
        iy_glb = iy_min + iy;
        for ( ix = 0; ix < nxl - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ix_glb = ix_min + ix;
            ii     = index( ix, iy );
            if ( !is_ghost_SR[ii] ) {
                ++ncells;

                // mark this node as being used
                if ( !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                // mark lowerright corner
                ii = index( ix + 1, iy );
                if ( !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                // mark upperleft corner
                ii = index( ix, iy + 1 );
                if ( !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                // mark upperright corner
                ii = index( ix + 1, iy + 1 );
                if ( !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
            }
        }
        nnodes_prev = nnodes;
    }
    nnodes += nb_extra_nodes;

    if( nparts == 1 ) {
        ATLAS_ASSERT( (periodic_x ? (nx+1)*ny_halo_NS : nx_halo_WE*ny_halo_NS )+ nb_extra_nodes == nnodes );
    }

    // define nodes and associated properties
    mesh.nodes().resize( nnodes );
    mesh::Nodes& nodes = mesh.nodes();
    auto ij            = array::make_view<idx_t, 2>(
        nodes.add( Field( "ij", array::make_datatype<idx_t>(), array::make_shape( nodes.size(), 2 ) ) ) );
    auto xy         = array::make_view<double, 2>( nodes.xy() );
    auto lonlat     = array::make_view<double, 2>( nodes.lonlat() );
    auto glb_idx    = array::make_view<gidx_t, 1>( nodes.global_index() );
    auto remote_idx = array::make_indexview<idx_t, 1>( nodes.remote_index() );
    auto part       = array::make_view<int, 1>( nodes.partition() );
    auto ghost      = array::make_view<int, 1>( nodes.ghost() );
    auto halo       = array::make_view<int, 1>( nodes.halo() );
    auto node_flags = array::make_view<int, 1>( nodes.flags() );
    auto lsm        = array::make_view<int, 1>(
        nodes.add( Field( "lsm", array::make_datatype<int>(), array::make_shape( nodes.size() ) ) ) );
    auto core = array::make_view<int, 1>(
        nodes.add( Field( "core", array::make_datatype<int>(), array::make_shape( nodes.size() ) ) ) );


    // define cells and associated properties
    mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
    if( include_south_pole ) {
        mesh.cells().add( new mesh::temporary::Triangle(), nxl-1 );
    }
    int quad_begin                                        = mesh.cells().elements( 0 ).begin();
    auto cells_part                                       = array::make_view<int, 1>( mesh.cells().partition() );
    auto cells_halo = array::make_view<int,1>( mesh.cells().halo());
    mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    auto flags = [&]( idx_t i ) { return Topology::view( node_flags( i ) ); };

    std::array<int,4> quad_nodes;
    int jcell = quad_begin;
    int inode, inode_nonghost, inode_ghost;

    // global indices for periodicity points
    std::vector<int> glb_idx_px;
    if ( periodic_x ) {
        glb_idx_px.resize( ny_halo_NS, -1 );
        inode = nx * ny_halo_NS;
        for ( iy = 0; iy < ny_halo_NS; iy++ ) {
            glb_idx_px[iy] = inode++;
        }
    }

    // loop over nodes and set properties
    ii             = 0;
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes

    int ix_pivot = nx / 2;
    bool patch = not orca.ghost( ix_pivot + 1, ny - 1 );
    for ( iy = 0; iy < nyl; iy++ ) {
        iy_glb = iy_min + iy;
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
        for ( ix = 0; ix < nxl; ix++ ) {
            ix_glb = ( ix_min + ix );  // don't take modulus here: periodicity points
                                       // have their own global index.
            // node properties
            if ( is_node_SR[ii] ) {
                // set node counter
                if ( is_ghost_SR[ii] ) {
                    inode = inode_ghost++;
                    ATLAS_ASSERT( inode < nnodes );
                }
                else {
                    inode = inode_nonghost++;
                    ATLAS_ASSERT( inode < nnodes_nonghost );
                }
                // global index
                if ( ix_glb < nx || not periodic_x ) {
                    ii_glb = local_idx_offset + iy_glb * local_idx_jstride + ix_glb;  // no periodic point
                }
                else {
                    ii_glb = glb_idx_px.at(iy_glb - iy_glb_min);
                }
                glb_idx( inode ) = ii_glb + 1;  // starting from 1

                // grid coordinates

                ij( inode, XX ) = ix_glb;
                ij( inode, YY ) = iy_glb;

                double _xy[2];

                // normal calculation
                orca.xy( ix_glb, iy_glb, _xy );
                normalise( _xy );

                xy( inode, LON ) = _xy[LON];
                xy( inode, LAT ) = _xy[LAT];

                // geographic coordinates by using projection
                lonlat( inode, LON ) = _xy[LON];
                lonlat( inode, LAT ) = _xy[LAT];

                // part
                part( inode ) = parts_SR[ii];
                // ghost nodes
                ghost( inode ) = is_ghost_SR[ii] || orca.ghost( ix_glb, iy_glb );

                if ( orca.ghost( ix_glb, iy_glb ) ) {
                    if ( not patch && iy_glb >= ny - 1 ) {  // ORCA2_T, ORCA025_T, ORCA12_T
                        int iy_pivot = ny - 1;
                        int ix_fold  = ix_pivot - ( ix_glb - ix_pivot );
                        int iy_fold  = iy_pivot - ( iy_glb - iy_pivot );
#if 0
                        // known to fail for ORCA1_F
                        // known to work for ORCA2_T
                        auto p_fold = orca.lonlat( ix_fold, iy_fold );
                        normalise( p_fold.data() );
                        auto p_this = orca.lonlat( ix_glb, iy_glb );
                        normalise( p_this.data() );
                        Log::info() << "[" << mpi::rank() << "]  (" << ix_glb << "," << iy_glb << ") " << p_this
                                    << "   (" << ix_fold << "," << iy_fold << ") " << p_fold << std::endl;
                        //   ATLAS_ASSERT( p_fold == p_this ); fails due to bugs in the actual data files ?
#endif
                        part( inode ) = partition( ix_fold, iy_fold );
                    }
                    else if ( patch && iy_glb >= ny ) {  // ORCA1_T
                        double iy_pivot = double( ny ) - 0.5;
                        int ix_fold     = 2 * ix_pivot - 1 - ix_glb;
                        int iy_fold     = int( 2. * iy_pivot ) - iy_glb;
#if 0
                        auto p_fold = orca.lonlat( ix_fold, iy_fold );
                        normalise( p_fold.data() );
                        auto p_this = orca.lonlat( ix_glb, iy_glb );
                        normalise( p_this.data() );
                        Log::info() << "[" << mpi::rank() << "] " << p_this << "   " << p_fold << std::endl;
                        ATLAS_ASSERT( p_fold == p_this );
#endif
                        part( inode ) = partition( ix_fold, iy_fold );
                    }
                    else if ( ix_glb >= nx ) {
                        part( inode ) = partition( nx - ix_glb, iy_glb );
                    }
                    else if ( ix_glb < 0 ) {
                        part( inode ) = partition( nx + ix_glb, iy_glb );
                    }
                    else {
                        // something here?
                    }
                }

                // flags
                flags( inode ).reset();
                if ( ghost( inode ) ) {
                    //                    ATLAS_DEBUG( "ghost: " << ix_glb << ", " << iy_glb << "  inode = " << inode
                    //                                           << " P = " << PointLonLat( _xy ) );
                    flags( inode ).set( Topology::GHOST );
                    remote_idx( inode ) = local_idx_SR[ii];
                    // change local index -- required for cells
                    local_idx_SR[ii] = inode;
                }
                else {
                    remote_idx( inode ) = -1;
                }

                if ( orca.land( ix_glb, iy_glb ) == 0 ) {
                    flags( inode ).set( Topology::LAND );
                }
                else {
                    flags( inode ).set( Topology::WATER );
                }
                lsm( inode )  = orca.water( ix_glb, iy_glb );
                core( inode ) = not orca.ghost( ix_glb, iy_glb );
                halo( inode ) = [&]() -> int {
                        if( ix_glb < 0 ) {
                            return -ix_glb;
                        }
                        else if( ix_glb > nx+1 ) {
                            return ix_glb - nx+1;
                        }
                        else if( iy_glb < 0 ) {
                            return -iy_glb;
                        }
                        else if( iy_glb >= ny ) {
                            int h = 1;
                            if ( patch && ix_glb < ix_pivot ) { // case of eg ORCA1_T
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
    int inode_south = ii;
    if( include_south_pole ) {
        inode = inode_south;
        PointXY p{orca.xy(0,0).x()+180.,-90.};
        lonlat(inode,LON) = p.x();
        lonlat(inode,LAT) = p.y();
        xy(inode,XX) = p.x();
        xy(inode,YY) = p.y();
        ij(inode,XX) = 0;
        ij(inode,YY) = iy_glb_min - 1;
        ghost(inode) = true;
        lsm(inode) = 0; // land
        flags(inode).set( Topology::LAND | Topology::SOUTH );
        halo(inode) = orca->haloSouth()+1;
        core(inode) = 0;
        remote_idx(inode) = -1;
        ++ii;
    }

    auto elem_flags_view = array::make_view<int, 1>( mesh.cells().flags() );
    auto elem_flags      = [&]( idx_t i ) { return Topology::view( elem_flags_view( i ) ); };

    // loop over nodes and define cells
    for ( iy = 0; iy < nyl - 1; iy++ ) {      // don't loop into ghost/periodicity row
        for ( ix = 0; ix < nxl - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ii = index( ix, iy );
            int ix_glb   = ix_min + ix;
            int iy_glb   = iy_min + iy;
            if ( !is_ghost_SR[ii] ) {
                    // define cell corners (local indices)
                    quad_nodes[0] = local_idx_SR[index( ix, iy )];          // lower left
                    quad_nodes[1] = local_idx_SR[index( ix + 1, iy )];      // lower right
                    quad_nodes[2] = local_idx_SR[index( ix + 1, iy + 1 )];  // upper right
                    quad_nodes[3] = local_idx_SR[index( ix, iy + 1 )];      // upper left

                    elem_flags( jcell ).reset();
                    if ( iy + 1 >= ny ) {
                        // std::swap( quad_nodes[2], quad_nodes[3] );
                    }

                    node_connectivity.set( jcell, quad_nodes.data() );
                    cells_part( jcell ) = part( quad_nodes[0] );

                    if ( iy_glb >= ny - 1 ) {
                        elem_flags( jcell ).set( Topology::GHOST );
                        if ( patch && ix_glb < ix_pivot ) {               // case of eg ORCA1_T
                            cells_part( jcell ) = part( quad_nodes[0] );  // lower left
                            elem_flags(jcell).unset( Topology::GHOST );
                        }
                        else {                                            // case of eg ORCA2_T
                            cells_part( jcell ) = part( quad_nodes[2] );  // upper right
                        }
                    }

                    bool elem_contains_water_point = [&] {
                        for ( idx_t inode: quad_nodes ) {
                            if ( flags( inode ).check( Topology::WATER ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    bool elem_contains_land_point = [&] {
                        for ( idx_t inode: quad_nodes ) {
                            if ( flags( inode ).check( Topology::LAND ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    cells_halo(jcell) = [&] {
                        int h=0;
                        for ( idx_t inode: quad_nodes ) {
                            h = std::max(h,halo(inode));
                        }
                        return h;
                    }();

                    if ( elem_contains_water_point ) {
                        elem_flags( jcell ).set( Topology::WATER );
                    }
                    if ( elem_contains_land_point ) {
                        elem_flags( jcell ).set( Topology::LAND );
                    }

                    ++jcell;
                }
        }
    }

    if( include_south_pole ) {
        int triag_begin  = mesh.cells().elements( 1 ).begin();
        int jcell = triag_begin;
        std::array<int,3> triag_nodes;
        for( int ix=0; ix<nxl-1; ++ix ) {
            ATLAS_ASSERT(jcell < mesh.cells().size());
            triag_nodes[0] = local_idx_SR[ index(ix,0) ];
            triag_nodes[2] = local_idx_SR[ index(ix+1,0) ];
            triag_nodes[1] = inode_south;
            node_connectivity.set(jcell,triag_nodes.data());
            elem_flags(jcell).set( Topology::LAND | Topology::GHOST );
            cells_halo(jcell) = [&] {
                int h=0;
                for( idx_t inode: triag_nodes ) {
                    h = std::max(h,halo(inode));
                }
                return h;
            }();
            cells_part(jcell) = 0;
            ++jcell;
        }
        ATLAS_ASSERT(jcell == mesh.cells().size());
    }

    generateGlobalElementNumbering( mesh );

    nodes.metadata().set( "parallel", true );
    nodes.metadata().set<size_t>( "nbRealPts", nnodes - nb_extra_nodes );
    nodes.metadata().set<size_t>( "NbVirtualPts", size_t( nb_extra_nodes ) );
}

OrcaMeshGenerator::OrcaMeshGenerator(const eckit::Parametrisation & config) {
    config.get("include_pole",include_pole_);
    config.get("force_include_south_pole",include_pole_);
    config.get("part", mypart_ = mpi::rank() );
    config.get("nb_parts", nparts_ = mpi::size() );
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
MeshGeneratorBuilder<OrcaMeshGenerator> __OrcaMeshGenerator( "orca" );
}

}  // namespace meshgenerator
}  // namespace atlas
