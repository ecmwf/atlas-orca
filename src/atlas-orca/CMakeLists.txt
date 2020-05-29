ecbuild_add_library(
    TARGET atlas-orca
    SOURCES Library.cc
    HEADER_DESTINATION "${INSTALL_INCLUDE_DIR}/atlas-orca"
    PUBLIC_LIBS eckit
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-orca>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)

