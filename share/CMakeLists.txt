file(GLOB orca_grids "atlas-orca/data/*.ascii")
install(
    FILES       ${orca_grids}
    DESTINATION share/atlas-orca/data
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/share/atlas-orca)
execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
    "${CMAKE_CURRENT_SOURCE_DIR}/atlas-orca/data"
    "${CMAKE_BINARY_DIR}/share/atlas-orca/data")

