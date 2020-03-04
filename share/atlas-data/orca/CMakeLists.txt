
set(files 
    "eorca025_F.ascii"
    "eorca025_T.ascii"
    "eorca025_U.ascii"
    "eorca025_V.ascii"
    "eorca1_F.ascii"
    "eorca1_T.ascii"
    "eorca1_U.ascii"
    "eorca1_V.ascii"
    "orca025_F.ascii"
    "orca025_T.ascii"
    "orca025_U.ascii"
    "orca025_V.ascii"
    "orca1_F.ascii"
    "orca1_T.ascii"
    "orca1_U.ascii"
    "orca1_V.ascii"
    "orca2_F.ascii"
    "orca2_T.ascii"
    "orca2_U.ascii"
    "orca2_V.ascii")

install(
    FILES ${files}
    DESTINATION "${CMAKE_BINARY_DIR}/share/${PROJECT_NAME}/orca"
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

foreach(_f ${files})
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${_f}" "${CMAKE_CURRENT_BINARY_DIR}/${_f}")
endforeach()

