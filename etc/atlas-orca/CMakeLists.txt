install(
    FILES       grid-specs.yaml
    DESTINATION etc/atlas-orca
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/etc/atlas-orca)
execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
    "${CMAKE_CURRENT_SOURCE_DIR}/grid-specs.yaml"
    "${CMAKE_BINARY_DIR}/etc/atlas-orca/grid-specs.yaml")

