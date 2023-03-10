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

add_custom_command( OUTPUT atlas-orca-cache.done
                    COMMAND ATLAS_ORCA_CACHING=1 ATLAS_CACHE_PATH=${CMAKE_BINARY_DIR}/share $<TARGET_FILE:atlas-orca-cache> --grid=all
                    COMMAND "${CMAKE_COMMAND}" -E touch atlas-orca-cache.done
                    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../../share/atlas/grids/orca/grids.yaml atlas-orca-cache
                    VERBATIM)

if ( ${ENABLE_ORCA_GRIDS_DOWNLOAD} )
  add_custom_target( atlas-orca-download ALL DEPENDS atlas-orca-cache.done )

  install( DIRECTORY ${CMAKE_BINARY_DIR}/share/atlas/grids/orca/v0/
           DESTINATION share/atlas/grids/orca/v0/
           DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
           FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ )
else()
  add_custom_target( atlas-orca-download DEPENDS atlas-orca-cache.done )

  install( DIRECTORY ${CMAKE_BINARY_DIR}/share/atlas/grids/orca/v0/
           DESTINATION share/atlas/grids/orca/v0/
           DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
           FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ )
endif()

add_subdirectory( atlas-orca-convert )
