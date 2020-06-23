
configure_file(version.h.in  version.h)
configure_file(version.cc.in version.cc)

ecbuild_add_library(
    TARGET atlas-orca
    SOURCES 
        Library.cc
        Library.h
        Orca.h
        Orca.cc
        OrcaGrid.h
        OrcaGrid.cc
        OrcaMeshGenerator.cc
        OrcaMeshGenerator.h
        ${CMAKE_CURRENT_BINARY_DIR}/version.cc
        ${CMAKE_CURRENT_BINARY_DIR}/version.h
    HEADER_DESTINATION "include/atlas-orca"
    PUBLIC_LIBS atlas
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-orca>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
       $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>)

