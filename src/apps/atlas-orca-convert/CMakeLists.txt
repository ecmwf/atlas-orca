# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

find_package( jedicmake REQUIRED )  # Prefer find modules from jedi-cmake
find_package( NetCDF COMPONENTS C CXX )

# Overcome Cray linking problem where the C++ library is not linked with the C library
list( APPEND NetCDF_CXX_EXTRA_LIBRARIES NetCDF::NetCDF_C )

ecbuild_add_option( FEATURE NETCDF
                    DESCRIPTION "Compile support to convert netcdf file to atlas-orca binary format"
                    REQUIRED_PACKAGES "NetCDF COMPONENTS C CXX" )

ecbuild_add_executable( TARGET atlas-orca-convert
                        SOURCES atlas-orca-convert.cc
                                AsciiReader.h
                                NetCDFReader.h
                                NetCDFReader.cc
                        LIBS atlas-orca )

target_compile_definitions( atlas-orca-convert PRIVATE ATLAS_ORCA_HAVE_NETCDF=${HAVE_NETCDF} )
if( HAVE_NETCDF )
    target_link_libraries( atlas-orca-convert NetCDF::NetCDF_CXX )
endif()
