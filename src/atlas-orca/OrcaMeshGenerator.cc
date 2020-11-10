
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

    int mypart  = options_.getInt( "part", mpi::rank() );
    int nparts  = options_.getInt( "nb_parts", mpi::size() );
    int nx      = rg.nx();
    int ny      = rg.ny();
    int ny_halo = ny + orca->haloNorth();

    bool periodic_x = true;  // options.get<bool>( "periodic_x" ) or rg.periodic();

// for asynchronous output
#if DEBUG_OUTPUT
    sleep( mypart );
#endif

    // this function should do the following:
    // - define nodes with
    //      mesh.nodes().resize(nnodes);
    //      mesh::Nodes& nodes = mesh.nodes();
    //    following properties should be defined:
    //      array::ArrayView<double,2> xy            ( nodes.xy() );
    //      array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
    //      array::ArrayView<int,   1> part          ( nodes.partition() );
    //      array::ArrayView<int,   1> ghost         ( nodes.ghost() );
    //      array::ArrayView<int,   1> flags         ( nodes.flags() );
    // - define cells (only quadrilaterals for now) with
    //      mesh.cells().add( new mesh::temporary::Quadrilateral(), nquads  );
    //    further define cells with
    //      array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index()
    //      );
    //      array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
    // - define connectivity with
    //      mesh::HybridElements::Connectivity& node_connectivity =
    //      mesh.cells().node_connectivity();
    //      node_connectivity.set( jcell, quad_nodes );
    //    where quad_nodes is a 4-element integer array containing the LOCAL
    //    indices of the nodes

    // Start with calculating number of quadrilaterals
    // The rule do determine if a cell belongs to a proc is the following: if the
    // lowerleft corner of the cell belongs to that proc.
    // so we loop over all gridpoints, select those that belong to the proc, and
    // determine the number of cells
    gidx_t ii_glb;  // global index
    int ncells;

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<idx_t> local_idx( nx * ny_halo, -1 );
    std::vector<int> current_idx( nparts, 0 );  // index counter for each proc

    auto partition = [&]( idx_t i, idx_t j ) -> int {
        j = std::min( ny - 1, j );
        i = std::min( nx - 1, i );
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
    for ( iy = 0; iy < ny_halo; iy++ ) {
        for ( ix = 0; ix < nx; ix++ ) {
            int p             = partition( ix, iy );
            local_idx[ii_glb] = current_idx[p]++;  // store local index on
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
    ix_max++;
    iy_max = std::min( ny_halo - 1, iy_max + 1 );
#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "SR = " << ix_min << ":" << ix_max << " x " << iy_min << ":" << iy_max << std::endl;
#endif

    // dimensions of surrounding rectangle (SR)
    int nxl = ix_max - ix_min + 1;
    int nyl = iy_max - iy_min + 1;

    // upper estimate for number of nodes
    nnodes_SR = nxl * nyl;
    ATLAS_DEBUG_VAR( nxl );
    ATLAS_DEBUG_VAR( nyl );
    ATLAS_DEBUG_VAR( ny );
    ATLAS_DEBUG_VAR( ny_halo );
    ATLAS_DEBUG_VAR( nnodes_SR );

    // partitions and local indices in SR
    std::vector<int> parts_SR( nnodes_SR, -1 );
    std::vector<int> local_idx_SR( nnodes_SR, -1 );
    std::vector<int> is_ghost_SR( nnodes_SR, 1 );
    ii = 0;  // index inside SR
    for ( iy = 0; iy < nyl; iy++ ) {
        iy_glb = iy_min + iy;  // global y-index
        for ( ix = 0; ix < nxl; ix++ ) {
            ix_glb       = ix_min + ix;  // global x-index
            parts_SR[ii] = partition( ix_glb, iy_glb );

            if ( ix_glb < nx && iy_glb < ny_halo ) {
                local_idx_SR[ii] = local_idx[iy_glb * nx + ix_glb];
                is_ghost_SR[ii]  = ( parts_SR[ii] != mypart );
            }
            else {
                local_idx_SR[ii] = -1;
                is_ghost_SR[ii]  = true;
                //                ATLAS_DEBUG( "ghost: " << ix_glb << "," << iy_glb );
            }
            ++ii;
        }
    }

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "parts_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << parts_SR[ii] << ",";
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "local_idx_SR = ";

    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << local_idx_SR[ii] << ",";
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "is_ghost_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << is_ghost_SR[ii] << ",";
    std::cout << std::endl;
#endif

    // vectors marking nodes that are necessary for this proc's cells
    std::vector<bool> is_node_SR( nnodes_SR, false );

    // determine number of cells and number of nodes
    nnodes     = 0;
    ncells     = 0;
    auto index = [&]( idx_t i, idx_t j ) { return j * nxl + i; };
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
    }
    ATLAS_DEBUG_VAR( nx * ny_halo );
    ATLAS_DEBUG_VAR( ( nx + 1 ) * ny_halo );
    ATLAS_DEBUG_VAR( nnodes );

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "nnodes = " << nnodes << std::endl;
    std::cout << "[" << mypart << "] : "
              << "is_node_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << is_node_SR[ii] << ",";
    std::cout << std::endl;
#endif

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
    int quad_begin                                        = mesh.cells().elements( 0 ).begin();
    auto cells_part                                       = array::make_view<int, 1>( mesh.cells().partition() );
    mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    auto flags = [&]( idx_t i ) { return Topology::view( node_flags( i ) ); };

    idx_t quad_nodes[4];
    int jcell = quad_begin;
    int inode, inode_nonghost, inode_ghost;

    // global indices for periodicity points
    inode = nx * ny_halo;
    std::vector<int> glb_idx_px( ny_halo, -1 );
    if ( periodic_x ) {
        for ( iy = 0; iy < ny_halo; iy++ ) {
            glb_idx_px[iy] = inode++;
        }
    }

    // loop over nodes and set properties
    ii             = 0;
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes

    bool patch = orca->core( nx / 2 + 1, ny - 1 );
    ATLAS_DEBUG_VAR( patch );

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
            //const double periodicity = ( ix_glb >= nx ) ? 360. : ( ix_glb < 0 ) ? -360. : 0.;
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
                if ( ix_glb < nx ) {
                    ii_glb = iy_glb * nx + ix_glb;  // no periodic point
                }
                else {
                    ii_glb = glb_idx_px[iy_glb];
                }
                glb_idx( inode ) = ii_glb + 1;  // starting from 1

                // grid coordinates

                ij( inode, XX ) = ix_glb;
                ij( inode, YY ) = iy_glb;

                double _xy[2];

                // normal calculation
                orca.xy( ix_glb, iy_glb, _xy );
                normalise( _xy );

                // testing Hack to get unique id (Don't ever do this)
                //                if ( !orca->core( ix_glb, iy_glb ) || iy_glb >= ny - 1 ) {
                //                    _xy[LAT] -= double( iy_glb ) * 0.00001;
                //                    _xy[LON] += double( ix_glb ) * 0.00001;
                //                }

                xy( inode, LON ) = _xy[LON];
                xy( inode, LAT ) = _xy[LAT];

                // geographic coordinates by using projection
                lonlat( inode, LON ) = _xy[LON];
                lonlat( inode, LAT ) = _xy[LAT];

                // part
                part( inode ) = parts_SR[ii];
                // ghost nodes
                ghost( inode ) = is_ghost_SR[ii] || !orca->core( ix_glb, iy_glb );

                if ( !orca->core( ix_glb, iy_glb ) ) {
                    int ix_pivot = nx / 2;
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

                if ( orca->lsm( ix_glb, iy_glb ) == 0 ) {
                    flags( inode ).set( Topology::LAND );
                }
                else {
                    flags( inode ).set( Topology::WATER );
                }
                lsm( inode )  = orca->lsm( ix_glb, iy_glb );
                core( inode ) = orca->core( ix_glb, iy_glb );
                halo( inode ) = 0;  // this is not yet correct
            }
            ++ii;
        }
    }
    auto elem_flags_view = array::make_view<int, 1>( mesh.cells().flags() );
    auto elem_flags      = [&]( idx_t i ) { return Topology::view( elem_flags_view( i ) ); };

    // loop over nodes and define cells
    for ( iy = 0; iy < nyl - 1; iy++ ) {      // don't loop into ghost/periodicity row
        for ( ix = 0; ix < nxl - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ii = index( ix, iy );
            if ( !is_ghost_SR[ii] ) {
                if ( ( ix_min + ix < nx - 1 || periodic_x ) && ( iy_min + iy < ny_halo - 1 ) ) {
                    // define cell corners (local indices)
                    quad_nodes[0] = local_idx_SR[index( ix, iy )];          // lower left
                    quad_nodes[1] = local_idx_SR[index( ix + 1, iy )];      // lower right
                    quad_nodes[2] = local_idx_SR[index( ix + 1, iy + 1 )];  // upper right
                    quad_nodes[3] = local_idx_SR[index( ix, iy + 1 )];      // upper left

                    elem_flags( jcell ).reset();
                    if ( iy + 1 >= ny ) {
                        // std::swap( quad_nodes[2], quad_nodes[3] );
                    }

                    node_connectivity.set( jcell, quad_nodes );
                    cells_part( jcell ) = part( quad_nodes[0] );


                    int ix_glb   = ix_min + ix;
                    int iy_glb   = iy_min + iy;
                    int ix_pivot = nx / 2;
                    if ( iy_glb >= ny - 1 ) {
                        elem_flags( jcell ).set( Topology::GHOST );
                        bool patch = orca->core( ix_pivot + 1, ny - 1 );
                        if ( patch && ix_glb < ix_pivot ) {               // case of eg ORCA1_T
                            cells_part( jcell ) = part( quad_nodes[0] );  // lower left
                        }
                        else {                                            // case of eg ORCA2_T
                            cells_part( jcell ) = part( quad_nodes[2] );  // upper right
                        }
                    }

                    bool elem_contains_water_point = [&] {
                        for ( idx_t jnode = 0; jnode < 4; ++jnode ) {
                            if ( flags( quad_nodes[jnode] ).check( Topology::WATER ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    bool elem_contains_land_point = [&] {
                        for ( idx_t jnode = 0; jnode < 4; ++jnode ) {
                            if ( flags( quad_nodes[jnode] ).check( Topology::LAND ) ) {
                                return true;
                            }
                        }
                        return false;
                    }();
                    if ( elem_contains_water_point ) {
                        elem_flags( jcell ).set( Topology::WATER );
                    }
                    if ( elem_contains_land_point ) {
                        elem_flags( jcell ).set( Topology::LAND );
                    }

//                    PointLonLat p0{lonlat( quad_nodes[0], LON ), lonlat( quad_nodes[0], LAT )};
//                    PointLonLat p1{lonlat( quad_nodes[1], LON ), lonlat( quad_nodes[1], LAT )};
//                    PointLonLat p2{lonlat( quad_nodes[2], LON ), lonlat( quad_nodes[2], LAT )};
//                    PointLonLat p3{lonlat( quad_nodes[3], LON ), lonlat( quad_nodes[3], LAT )};
//                    double two_degrees = geometry.distance( PointLonLat{0, 0}, PointLonLat{3, 0} );
//                    double d01         = geometry.distance( p0, p1 );
//                    double d02         = geometry.distance( p0, p2 );
//                    double d03         = geometry.distance( p0, p3 );
//                    double d12         = geometry.distance( p1, p2 );
//                    double d13         = geometry.distance( p1, p3 );
//                    double d23         = geometry.distance( p2, p3 );
//                    double d           = std::max( {d01, d02, d03, d12, d13, d23} );
//                    static int count   = 1;
//                    if ( count == 1 ) {
//                        ATLAS_DEBUG_VAR( two_degrees );
//                        count++;
//                    }
//                    if ( d > two_degrees ) {
//                        ATLAS_DEBUG( "invalid element " << PointLonLat( ix_glb, iy_glb ) << "  " << d );

//                        Topology::set( elem_flags( jcell ), Topology::GHOST );
//                    }
#if DEBUG_OUTPUT_DETAIL
                    std::cout << "[" << mypart << "] : "
                              << "New quad " << 1 + jcell << "\t";
                    std::cout << ": " << 1 + quad_nodes[0] << "," << 1 + quad_nodes[1] << "," << 1 + quad_nodes[2]
                              << "," << 1 + quad_nodes[3] << std::endl;
#endif

                    ++jcell;
                }
            }
        }
    }

#if DEBUG_OUTPUT
    // list nodes
    for ( inode = 0; inode < nnodes; inode++ ) {
        std::cout << "[" << mypart << "] : "
                  << " node " << inode << ": ghost = " << ghost( inode ) << ", glb_idx = " << glb_idx( inode ) - 1
                  << ", part = " << part( inode ) << ", lon = " << lonlat( inode, 0 )
                  << ", lat = " << lonlat( inode, 1 ) << ", remote_idx = " << remote_idx( inode ) << std::endl;
    }

    int* cell_nodes;
    for ( jcell = 0; jcell < ncells; jcell++ ) {
        std::cout << "[" << mypart << "] : "
                  << " cell " << jcell << ": " << node_connectivity( jcell, 0 ) << "," << node_connectivity( jcell, 1 )
                  << "," << node_connectivity( jcell, 2 ) << "," << node_connectivity( jcell, 3 ) << std::endl;
    }
#endif

    generateGlobalElementNumbering( mesh );

    nodes.metadata().set( "parallel", true );
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
