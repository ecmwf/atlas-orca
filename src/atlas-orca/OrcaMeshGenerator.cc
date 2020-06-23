
#include "OrcaMeshGenerator.h"
#include "OrcaGrid.h"

#include <tuple>

#include "atlas/array/MakeView.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/StructuredMeshGenerator.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/NormaliseLongitude.h"

#include "eckit/utils/Hash.h"

namespace atlas {
namespace meshgenerator {

void OrcaMeshGenerator::generate( const Grid& irreg_grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    util::Config config;
    //config.set( "patch_pole", false );
    StructuredMeshGenerator regular_mesh{config};

    OrcaGrid orca{irreg_grid};
    ATLAS_ASSERT( orca );

    // TODO: This hard-coding should instead go to the data-file, and become queried via the OrcaGrid
    std::map<std::string, bool> patch;
    std::vector<std::tuple<std::string, bool>> resolutions = {{"025", false}, {"12", false}, {"1", true}, {"2", false}};
    std::vector<std::string> prefixes                      = {"orca", "eorca"};
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
    StructuredGrid::YSpace yspace{grid::LinearSpacing{{-80., 90.}, orca.ny(), not patch[orca.name()]}};
    // Periodic xspace
    StructuredGrid::XSpace xspace{grid::LinearSpacing{{0., 360.}, orca.nx(), false}};
    StructuredGrid reg_grid{xspace, yspace};

    regular_mesh.generate( reg_grid, distribution, mesh );

    reassign_coordinates( orca, mesh );
}

void OrcaMeshGenerator::generate( const Grid& irreg_grid, Mesh& mesh ) const {
    grid::Partitioner partitioner( "serial" );
    grid::Distribution distribution{irreg_grid, partitioner};
    generate( irreg_grid, distribution, mesh );
}

namespace {
std::tuple<double, double> fetch_ghost( const OrcaGrid& orca, const gidx_t glid ) {
    idx_t i       = 0;
    idx_t j       = glid - orca.size();
    PointLonLat p = orca.lonlat( i, j );

    return {p.lon() + 360, p.lat()};
}

std::tuple<double, double> fetch_inner( const OrcaGrid& orca, const gidx_t glid ) {
    idx_t i       = glid % orca.nx();
    idx_t j       = glid / orca.nx();
    PointLonLat p = orca.lonlat( i, j );

    if ( i >= orca.nx() / 2 ) {
        util::NormaliseLongitude normalise{orca.lonlat( 0, j ).lon() + 360. / orca.nx()};
        return {normalise( p.lon() ), p.lat()};
    }
    else {
        util::NormaliseLongitude normalise{orca.lonlat( 0, j ).lon()};
        return {normalise( p.lon() ), p.lat()};
    }
}

}  // namespace

void OrcaMeshGenerator::reassign_coordinates( const OrcaGrid& orca, Mesh& mesh ) const {
    auto lonlat  = array::make_view<double, 2>( mesh.nodes().lonlat() );
    auto glb_idx = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );

    for ( idx_t n = 0; n < mesh.nodes().size(); ++n ) {
        gidx_t g = glb_idx[n] - 1;

        std::tie( lonlat( n, 0 ), lonlat( n, 1 ) ) =
            ( orca.size() <= g ) ? fetch_ghost( orca, g ) : fetch_inner( orca, g );

        ATLAS_ASSERT( g < mesh.nodes().size() );
    }
}

void OrcaMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "OrcaMeshGenerator" );
}

namespace {
MeshGeneratorBuilder<OrcaMeshGenerator> __OrcaMeshGenerator( "orca" );
}

}  // namespace meshgenerator
}  // namespace atlas
