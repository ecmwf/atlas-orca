configure_file(version.h.in  version.h)
configure_file(version.cc.in version.cc)

ecbuild_add_library(
    TARGET atlas-orca
    SOURCES 
        Library.cc
        Library.h
        grid/Orca.h
        grid/Orca.cc
        grid/OrcaGrid.h
        grid/OrcaGrid.cc
        meshgenerator/OrcaMeshGenerator.cc
        meshgenerator/OrcaMeshGenerator.h
        meshgenerator/FixupMesh.cc
        meshgenerator/FixupMesh.h
        meshgenerator/fixups/FixupMesh_ORCA.cc
        meshgenerator/fixups/FixupMesh_ORCA.h
        meshgenerator/fixups/FixupMesh_eORCA1.cc
        meshgenerator/fixups/FixupMesh_eORCA1.h
        meshgenerator/fixups/FixupMesh_ORCA2.cc
        meshgenerator/fixups/FixupMesh_ORCA2.h
        ${CMAKE_CURRENT_BINARY_DIR}/version.cc
        ${CMAKE_CURRENT_BINARY_DIR}/version.h
    HEADER_DESTINATION "include/atlas-orca"
    PUBLIC_LIBS atlas
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-orca>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
       $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>)

