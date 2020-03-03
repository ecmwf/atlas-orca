ecbuild_add_library(
    TARGET atlas-data
    SOURCES Library.cc
    HEADER_DESTINATION "${INSTALL_INCLUDE_DIR}/atlas-data"
    PUBLIC_LIBS eckit
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-data>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)

