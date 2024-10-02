/*
 * (C) Copyright 2021- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <bitset>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "eckit/filesystem/PathName.h"
#include "eckit/log/Bytes.h"
#include "eckit/types/Fraction.h"

#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"

#include "AsciiReader.h"
#include "NetCDFReader.h"
#include "atlas-orca/util/AtlasIOReader.h"
#include "atlas-orca/util/OrcaData.h"


namespace atlas::orca {

//----------------------------------------------------------------------------------------------------------------------

struct Tool : public atlas::AtlasTool {
    bool serial() override { return true; }
    int execute( const Args& args ) override;
    std::string briefDescription() override { return "Create binary grid data files "; }
    std::string usage() override {
        return name() + " <file> --name=NAME --arrangement={F,T,U,V,W} [OPTION]... [--help,-h]";
    }
    std::string longDescription() override {
        return "Create binary grid data files \n"
               "\n"
               "       <file>: input file";
    }

    Tool( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "name", "Output grid name" ) );
        add_option( new SimpleOption<std::string>( "arrangement", "Arrangement in Arakawa C-grid: F,T,U,V,W" ) );
        add_option( new SimpleOption<std::string>( "input-format",
                                                   "netcdf, ascii-v1, ascii-v2, atlas-io; default: ascii-v2" ) );
        add_option( new SimpleOption<std::string>( "compression",
                                                   "Data compression: none, lz4, aec, ... (see eckit support)'" ) );
        add_option(
            new SimpleOption<std::string>( "output", "Output file path; default: <name>_<arrangement>.atlas" ) );
        add_option( new Separator( "Advanced" ) );
        add_option( new SimpleOption<bool>( "verbose", "Print verbose output" ) );
        add_option( new SimpleOption<double>(
            "diagonal-factor", "Diagonal factor used to limit valid elements, with respect to nominal resolution" ) );
        add_option( new VectorOption<double>(
            "pivot", "Pivot point of North fold, if cannot automatically be determined", 2, "," ) );
        add_option( new SimpleOption<bool>( "yaml", "Output spec instead of data" ) );
    }
};

//------------------------------------------------------------------------------------------------------


int Tool::execute( const Args& args ) {
    // User sanity checks
    if ( args.count() == 0 ) {
        Log::error() << "No file specified." << std::endl;
        help( std::cout );
        return failed();
    }
    if ( args.count() > 1 ) {
        Log::error() << "Only one file can be specified." << std::endl;
        help( std::cout );
        return failed();
    }
    std::string input{ args( 0 ) };
    if ( input.find( "http" ) != 0 ) {
        eckit::PathName file( input );
        if ( !file.exists() ) {
            Log::error() << "File does not exist: " << file << std::endl;
            return failed();
        }
    }

    bool yaml_output = args.getBool( "yaml", false );

    std::string name       = args.getString( "name", "unnamed" );
    std::string P          = args.getString( "arrangement", "P" );
    std::string outputfile = args.getString( "output", name + "_" + P + ".atlas" );

    std::string input_format;
    if ( not args.get( "input-format", input_format ) ) {
        std::string extension = input.substr( input.find_last_of( '.' ) );
        if ( extension == ".nc" ) {
            input_format = "netcdf";
        }
        else if ( extension == ".ascii" ) {
            input_format = "ascii-v2";
        }
        else if ( extension == ".atlas" ) {
            input_format = "atlas-io";
        }
        else {
            Log::warning() << "Could not determine input-format from extension " << extension << std::endl;
            return failed();
        }
    }

    if ( yaml_output ) {
        Log::info().reset();
    }

    OrcaData data;

    if ( input_format.find( "netcdf" ) == 0 ) {
        NetCDFReader{ args }.read( input, data );
    }
    else if ( input_format.find( "ascii-v" ) == 0 ) {
        std::string version_str = input_format.substr( 7 );
        int version             = std::stoi( version_str );
        util::Config reader_config( args );
        reader_config.set( "version", version );
        AsciiReader reader( reader_config );
        reader.read( input, data );
    }
    else if ( input_format == "atlas-io" ) {
        AtlasIOReader{ util::NoConfig() }.read( input, data );
    }
    else {
        Log::warning() << "Unknown input_format \"" << input_format << "\"" << std::endl;
        Log::warning() << "Supported: atlas-io, ascii-v1, ascii-v2" << std::endl;
        return failed();
    }

    std::string uid = data.computeUid( args );

    if ( not yaml_output ) {
        auto invalid_element_statistics = data.detectInvalidElements( args );
        // Print summary
        const int halo_N = data.halo[HALO_NORTH];
        const int halo_W = data.halo[HALO_WEST];
        const int halo_S = data.halo[HALO_SOUTH];
        const int halo_E = data.halo[HALO_EAST];
        int nx           = data.dimensions[0] - halo_W - halo_E;

        Log::info() << "resolution       : " << eckit::Fraction{ 360, nx } << " degrees" << std::endl;
        Log::info() << "dimensions       : "
                    << "[" << data.dimensions[0] << "," << data.dimensions[1] << "]" << std::endl;
        Log::info() << "halo [N,W,S,E]   : "
                    << "[" << halo_N << "," << halo_W << "," << halo_S << "," << halo_E << "]" << std::endl;
        Log::info() << "invalid elements : " << invalid_element_statistics.invalid_elements << std::endl;
        if ( invalid_element_statistics.invalid_elements > 0 ) {
            Log::info() << "    quad2d       : " << invalid_element_statistics.invalid_quads_2d << std::endl;
            Log::info() << "    quad3d       : " << invalid_element_statistics.invalid_quads_3d << std::endl;
            Log::info() << "    large diag   : " << invalid_element_statistics.diagonal_too_large << std::endl;
        }
        Log::info() << "uid              : " << uid << std::endl;
        Log::info() << "pivot            : " << data.pivot << std::endl;
    }

    if ( not yaml_output ) {
        ATLAS_TRACE( "Write data file" );
        auto length = data.write( outputfile, args );
        Log::info() << "Written " << eckit::Bytes( static_cast<double>( length ) ) << " to file " << outputfile
                    << std::endl;
    }

    if ( yaml_output ) {
        std::string P     = args.getString( "arrangement", "P" );
        std::ostream& out = std::cout;
        out << name << "_" << P << ": &" << name << "_" << P << std::endl;
        out << "    type: ORCA" << std::endl;
        out << "    orca_arrangement: " << P << std::endl;
        out << "    orca_name: " << name << std::endl;
        out << "    dimensions: " << data.dimensions << std::endl;
        out << "    uid: " << uid << std::endl;
        out << "    data: {{location}}/" << outputfile << std::endl;
        out << uid << ": *" << name << "_" << P << std::endl;
    }

    return success();
}

}  // namespace atlas::orca


//------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    atlas::orca::Tool tool( argc, argv );
    return tool.start();
}
