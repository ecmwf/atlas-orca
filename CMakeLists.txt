# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

find_package(ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project(atlas-orca LANGUAGES CXX)

find_package(eckit 1.16.1 REQUIRED)
find_package(atlas 0.25.0 REQUIRED)

if( NOT eckit_HAVE_LZ4 )
    ecbuild_warn( "eckit is required to have LZ4 support to decompress atlas-orca binary data files. "
                  "Continuing regardless with reduced tests" )
endif()

ecbuild_add_option( FEATURE ORCA_GRIDS_DOWNLOAD
                    DEFAULT OFF
                    DESCRIPTION "Download ORCA grids during build" )

atlas_create_plugin( atlas-orca
    URL       https://sites.ecmwf.int/docs/atlas
    NAMESPACE int.ecmwf )

add_subdirectory(etc)
add_subdirectory(share)
add_subdirectory(src)

ecbuild_install_project(NAME ${PROJECT_NAME})

