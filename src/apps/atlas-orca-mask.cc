/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/log/JSON.h"

#include "atlas/io/atlas-io.h"

#if ATLAS_ORCA_HAVE_ECKIT_CODEC
#include "eckit/codec/codec.h"
#else
namespace eckit::codec {
using RecordReader = atlas::io::RecordReader;
}
#endif

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/SpecRegistry.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/AtlasTool.h"

#include "atlas-orca/util/Enums.h"
#include "atlas-orca/util/Flag.h"
#include "atlas-orca/util/OrcaDataFile.h"

namespace atlas::orca {

namespace {

constexpr size_t MAX_REPORTED_INVALID_VALUES = 20;
constexpr char MASK_FORMAT[]                  = "atlas-mask";
constexpr int MASK_VERSION                    = 1;

struct GridData {
    Grid::Spec grid_spec;
    bool has_grid_spec{ false };
    std::array<std::int32_t, 2> dimensions;
    std::array<std::int32_t, 4> halo;
    std::vector<std::byte> flags;
};

Grid::Spec resolve_grid_spec( const std::string& specification ) {
    const auto first = specification.find_first_not_of( " \t\r\n" );
    if ( first != std::string::npos && specification[first] == '{' ) {
        std::istringstream stream( specification );
        return util::Config( stream );
    }

    if ( grid::SpecRegistry::has( specification ) ) {
        return grid::SpecRegistry::get( specification );
    }

    std::string identifier = specification;
    std::transform( identifier.begin(), identifier.end(), identifier.begin(),
                    []( unsigned char character ) { return static_cast<char>( std::tolower( character ) ); } );
    if ( grid::SpecRegistry::has( identifier ) ) {
        return grid::SpecRegistry::get( identifier );
    }

    std::transform( identifier.begin(), identifier.end(), identifier.begin(),
                    []( unsigned char character ) { return static_cast<char>( std::toupper( character ) ); } );
    if ( !identifier.empty() && identifier.front() == 'E' ) {
        identifier.front() = 'e';
    }
    if ( grid::SpecRegistry::has( identifier ) ) {
        return grid::SpecRegistry::get( identifier );
    }

    throw std::invalid_argument( "Unknown grid name or UID '" + specification + "'" );
}

std::string grid_name( const Grid::Spec& spec ) {
    const std::string name        = spec.getString( "orca_name", "" );
    const std::string arrangement = spec.getString( "orca_arrangement", "" );
    return name.empty() || arrangement.empty() ? spec.getString( "uid", "grid" ) : name + "_" + arrangement;
}

size_t mismatch_limit( const std::string& value ) {
    if ( value == "unlimited" ) {
        return std::numeric_limits<size_t>::max();
    }
    if ( value.empty() || value.front() == '-' ) {
        throw std::invalid_argument( "Invalid mismatch limit" );
    }

    size_t parsed     = 0;
    const auto limit = std::stoull( value, &parsed );
    if ( parsed != value.size() || limit > std::numeric_limits<size_t>::max() ) {
        throw std::invalid_argument( "Invalid mismatch limit" );
    }
    return static_cast<size_t>( limit );
}

bool is_grid_halo( size_t index, const std::array<std::int32_t, 2>& dimensions,
                   const std::array<std::int32_t, 4>& halo ) {
    const size_t i = index % dimensions[0];
    const size_t j = index / dimensions[0];
    return i < halo[HALO_WEST] || i >= dimensions[0] - halo[HALO_EAST] || j < halo[HALO_SOUTH] ||
           j >= dimensions[1] - halo[HALO_NORTH];
}

void read_grid_data( const std::string& source, GridData& data ) {
    eckit::codec::RecordReader reader( OrcaDataFile{ source } );
    reader.read( "dimensions", data.dimensions );
    reader.read( "halo", data.halo );
    reader.read( "flags", data.flags );
    reader.wait();
}

void write_mask( const std::string& path, const GridData& data ) {
    std::vector<int> mask;
    mask.reserve( data.flags.size() );
    for ( const auto& flag : data.flags ) {
        mask.emplace_back( Flag{ flag }.test( Flag::WATER ) );
    }

    eckit::LocalConfiguration compression;
    compression.set( "compression", "lz4" );

    io::RecordWriter writer;
    writer.set( "format", MASK_FORMAT );
    writer.set( "version", MASK_VERSION );
    writer.set( "grid", data.grid_spec.json( eckit::JSON::Formatting::compact() ) );
    writer.set( "mask", io::ref( mask ), compression );
    writer.write( path, io::Mode::write );
}

bool read_mask( const std::string& path, Grid::Spec& grid_spec, std::vector<int>& mask ) {
    io::InputFileStream stream( path );
    io::Record record;
    record.read( stream );

    const bool has_format  = record.has( "format" );
    const bool has_version = record.has( "version" );
    if ( has_format != has_version ) {
        Log::error() << "Mask record '" << path << "' must contain both 'format' and 'version'." << std::endl;
        return false;
    }

    std::string format;
    int version = 0;
    eckit::codec::RecordReader reader( path );
    if ( has_format ) {
        reader.read( "format", format );
        reader.read( "version", version );
    }
    std::string encoded_grid;
    reader.read( "grid", encoded_grid );
    reader.read( "mask", mask );
    reader.wait();

    if ( has_format ) {
        if ( format != MASK_FORMAT ) {
            Log::error() << "Mask record '" << path << "' has unsupported format '" << format << "'." << std::endl;
            return false;
        }
        if ( version != MASK_VERSION ) {
            Log::error() << "Mask record '" << path << "' has unsupported " << MASK_FORMAT << " version " << version
                         << "." << std::endl;
            return false;
        }

        std::istringstream stream( encoded_grid );
        grid_spec = util::Config( stream );
    }
    else if ( !encoded_grid.empty() ) {
        grid_spec = resolve_grid_spec( encoded_grid );
    }

    if ( encoded_grid.empty() ) {
        Log::error() << "Mask record '" << path << "' contains an empty grid name." << std::endl;
        return false;
    }
    return true;
}

void write_gmsh( const std::string& path, const Grid::Spec& grid_spec, const std::vector<int>& mask_diff ) {
    auto mesh  = Mesh( Grid( grid_spec ) );
    auto fs    = functionspace::NodeColumns( mesh );
    auto field = fs.createField<int>( option::name( "mask_diff" ) );
    auto view  = array::make_view<int, 1>( field );

    ATLAS_ASSERT_MSG( static_cast<size_t>( view.size() ) == mask_diff.size(),
                      "Mesh node count does not match mask size" );
    for ( idx_t index = 0; index < view.size(); ++index ) {
        view( index ) = mask_diff[index];
    }

    output::Gmsh gmsh( path, util::Config( "coordinates", "xyz" )( "ghost", true ) );
    gmsh.write( mesh );
    gmsh.write( field );
}

}  // namespace

class MaskTool : public atlas::AtlasTool {
public:
    MaskTool( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "grid", "ORCA grid name, UID, or JSON specification" ) );
        add_option( new SimpleOption<std::string>( "data", "ORCA .atlas data file or URL" ) );
        add_option( new SimpleOption<std::string>( "mask", "Mask file to validate" ) );
        add_option( new SimpleOption<std::string>( "output", "Output mask or Gmsh path" ) );
        add_option( new SimpleOption<bool>( "output-gmsh", "Write Gmsh differences when validation fails" ) );
        add_option(
            new SimpleOption<std::string>( "filter", "Reported mismatch ownership: all, halo (ghost), or owned" ) );
        add_option( new SimpleOption<std::string>( "limit", "Maximum reported mismatches, or unlimited" ) );
    }

    bool serial() override { return true; }
    std::string briefDescription() override { return "Create and validate ORCA water masks"; }
    std::string usage() override {
                return name() + " create --grid=GRID [--data=PATH] [--output=PATH]\n"
                         + name() + " validate (--grid=GRID [--data=PATH] | --data=PATH) --mask=PATH [--filter=all|halo|owned] "
                "[--limit=N|unlimited] [--output-gmsh] [--output=PATH] [--help,-h]";
    }
    std::string longDescription() override {
        return "Create a standalone integer water mask from a registered ORCA grid or .atlas data file, or validate "
                             "an existing mask against one. GRID may be a registered name, UID, or JSON specification; --data "
                             "overrides its data source.";
    }

    int execute( const Args& args ) override;

private:
    bool loadGridData( const Args& args, bool require_grid, GridData& data );
    int createMask( const Args& args, const GridData& data );
    int validateMask( const Args& args, const GridData& data );
};

bool MaskTool::loadGridData( const Args& args, bool require_grid, GridData& data ) {
    std::string grid_specification;
    std::string data_path;
    const bool has_grid = args.get( "grid", grid_specification );
    const bool has_data = args.get( "data", data_path );
    if ( !has_grid && !has_data ) {
        Log::error() << "Specify --grid, --data, or both." << std::endl;
        return false;
    }
    if ( require_grid && !has_grid ) {
        Log::error() << "--grid is required when creating a mask." << std::endl;
        return false;
    }

    if ( has_grid ) {
        ATLAS_DEBUG_VAR(grid_specification);
        data.grid_spec     = resolve_grid_spec( grid_specification );
        data.has_grid_spec = true;
        if ( !has_data ) {
            ATLAS_DEBUG_VAR(data.grid_spec);
            if ( !data.grid_spec.has( "data" ) ) {
                Log::error() << "Grid '" << grid_name( data.grid_spec ) << "' has no data source; specify --data."
                             << std::endl;
                return false;
            }
            data_path = data.grid_spec.getString( "data" );
            ATLAS_DEBUG_VAR(data_path);
        }
    }

    read_grid_data( data_path, data );
    return true;
}

int MaskTool::createMask( const Args& args, const GridData& data ) {
    const std::string output_path = args.getString( "output", "mask_" + grid_name( data.grid_spec ) + ".atlas" );
    write_mask( output_path, data );
    Log::info() << "Wrote " << data.flags.size() << " water-mask values to '" << output_path << "'." << std::endl;
    return success();
}

int MaskTool::validateMask( const Args& args, const GridData& data ) {
    std::string mask_path;
    if ( !args.get( "mask", mask_path ) || mask_path.empty() ) {
        Log::error() << "--mask is required for validation." << std::endl;
        return failed();
    }

    const std::string filter = args.getString( "filter", "all" );
    if ( filter != "all" && filter != "halo" && filter != "owned" ) {
        Log::error() << "Invalid filter '" << filter << "'; expected all, halo, or owned." << std::endl;
        return failed();
    }

    size_t report_limit;
    try {
        report_limit = mismatch_limit( args.getString( "limit", "20" ) );
    }
    catch ( const std::exception& ) {
        Log::error() << "Invalid limit; expected a non-negative integer or unlimited." << std::endl;
        return failed();
    }

    Grid::Spec mask_grid_spec;
    std::vector<int> mask;
    if ( !read_mask( mask_path, mask_grid_spec, mask ) ) {
        return failed();
    }

    if ( data.flags.size() != mask.size() ) {
        Log::error() << "Array size mismatch: grid flags contain " << data.flags.size() << " values, but mask contains "
                     << mask.size() << " values." << std::endl;
        return failed();
    }

    std::vector<int> mask_diff( data.flags.size() );
    size_t mismatches       = 0;
    size_t ghost_mismatches = 0;
    size_t halo_mismatches  = 0;
    size_t reported         = 0;
    size_t invalid_values   = 0;
    for ( size_t index = 0; index < data.flags.size(); ++index ) {
        const bool grid_water = Flag{ data.flags[index] }.test( Flag::WATER );
        if ( mask[index] != 0 && mask[index] != 1 ) {
            mask_diff[index] = -1;
            if ( invalid_values < MAX_REPORTED_INVALID_VALUES ) {
                Log::error() << "Invalid mask value at index " << index << ": " << mask[index]
                             << " (expected 0 or 1)" << std::endl;
            }
            ++invalid_values;
        }
        else if ( grid_water != static_cast<bool>( mask[index] ) ) {
            const bool ghost        = Flag{ data.flags[index] }.test( Flag::GHOST );
            const bool in_grid_halo = is_grid_halo( index, data.dimensions, data.halo );
            const bool report = filter == "all" || ( filter == "halo" && ghost ) || ( filter == "owned" && !ghost );
            mask_diff[index] = grid_water ? 3 : 2;
            if ( report && reported < report_limit ) {
                Log::error() << "Mismatch at index " << index << ": grid water=" << grid_water
                             << ", mask=" << mask[index] << ", i=" << std::setw( 4 ) << index % data.dimensions[0]
                             << ", j=" << std::setw( 4 ) << index / data.dimensions[0]
                             << ", ghost=" << ( ghost ? "true" : "false" )
                             << ", grid-halo=" << ( in_grid_halo ? "true" : "false" ) << std::endl;
                ++reported;
            }
            ++mismatches;
            ghost_mismatches += ghost;
            halo_mismatches += in_grid_halo;
        }
        else {
            mask_diff[index] = grid_water ? 1 : 0;
        }
    }

    if ( invalid_values != 0 || mismatches != 0 ) {
        Log::error() << "Validation failed: " << mismatches << " mismatches and " << invalid_values
                     << " invalid mask values among " << data.flags.size() << " entries; " << ghost_mismatches
                     << " mismatches are ghost points, of which " << halo_mismatches
                     << " are in the rectangular grid halo." << std::endl;
        if ( args.getBool( "output-gmsh", false ) ) {
            const Grid::Spec& grid_spec = data.has_grid_spec ? data.grid_spec : mask_grid_spec;
            const std::string output_path = args.getString( "output", "mask_diff.msh" );
            Log::info() << "Writing Gmsh mask differences to '" << output_path << "'." << std::endl;
            write_gmsh( output_path, grid_spec, mask_diff );
        }
        return failed();
    }

    Log::info() << "Validation passed: all " << data.flags.size() << " water-mask values are identical." << std::endl;
    return success();
}

int MaskTool::execute( const Args& args ) {
    if ( args.count() != 1 ) {
        Log::error() << "Expected a command: create or validate." << std::endl;
        help( std::cout );
        return failed();
    }

    const std::string command = args( 0 );
    if ( command != "create" && command != "validate" ) {
        Log::error() << "Unknown command '" << command << "'; expected create or validate." << std::endl;
        return failed();
    }

    GridData data;
    if ( !loadGridData( args, command == "create", data ) ) {
        return failed();
    }
    return command == "create" ? createMask( args, data ) : validateMask( args, data );
}

}  // namespace atlas::orca

int main( int argc, char** argv ) {
    atlas::orca::MaskTool tool( argc, argv );
    return tool.start();
}