# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


set( files       grids.yaml )
set( destination share/atlas/grids/orca )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination})

foreach( file ${files} )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_if_different"
        "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
        "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()


# Retrieve data files
# -------------------

add_custom_command( OUTPUT atlas-orca-retrieve-orca-data.done
                    COMMAND ATLAS_ORCA_CACHING=1 ATLAS_CACHE_PATH=${CMAKE_BINARY_DIR}/share $<TARGET_FILE:atlas-orca-cache> --grid=all
                    COMMAND "${CMAKE_COMMAND}" -E touch atlas-orca-retrieve-orca-data.done
                    DEPENDS ${PROJECT_SOURCE_DIR}/share/atlas/grids/orca/grids.yaml atlas-orca-cache
                    VERBATIM
                    COMMENT "Retrieve atlas-orca grid data for all grids present in ${PROJECT_SOURCE_DIR}/share/atlas/grids/orca/grids.yaml" )

if ( HAVE_RETRIEVE_ORCA_DATA )
  add_custom_target( atlas-orca-retrieve-orca-data ALL DEPENDS atlas-orca-retrieve-orca-data.done )
else()
  add_custom_target( atlas-orca-retrieve-orca-data DEPENDS atlas-orca-retrieve-orca-data.done )
endif()


# Installation
# ------------

if ( HAVE_INSTALL_ORCA_DATA )
    install(
        DIRECTORY ${CMAKE_BINARY_DIR}/share/atlas/grids/orca/
        DESTINATION ${destination}
        DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
        FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ )
else()
    install(
        FILES       ${files}
        DESTINATION ${destination}
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()
