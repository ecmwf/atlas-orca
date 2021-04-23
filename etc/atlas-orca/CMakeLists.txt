set( files       grids.yaml config.yaml )
set( destination etc/atlas-orca )

install(
    FILES       ${files}
    DESTINATION ${destination}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination})

foreach( file ${files} )
  execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
      "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
      "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()
