
configure_file(atlas_orca_config.h.in atlas_orca_config.h)
configure_file(atlas_orca_config.cc.in atlas_orca_config.cc)

ecbuild_add_library(
    TARGET atlas_orca
    SOURCES Library.cc ${CMAKE_CURRENT_BINARY_DIR}/atlas_orca_config.cc
    HEADER_DESTINATION "${INSTALL_INCLUDE_DIR}/atlas-orca"
    PUBLIC_LIBS eckit
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-orca>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
       $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>)

