# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable( TARGET atlas-orca-cache
                        SOURCES atlas-orca-cache.cc
                        LIBS atlas-orca )

add_subdirectory( atlas-orca-convert )

set( atlas_orca_scripts
    atlas-orca-recreate-grids
    atlas-orca-upload-grids )

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
foreach( script ${atlas_orca_scripts} )
  execute_process( COMMAND "${CMAKE_COMMAND}" -E create_symlink
                "${CMAKE_CURRENT_SOURCE_DIR}/${script}"
                "${CMAKE_BINARY_DIR}/bin/${script}" )
endforeach()
