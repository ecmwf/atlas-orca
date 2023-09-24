# (C) Copyright 2021- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

configure_file(version.cc.in version.cc)

if ( HAVE_ORCA_IN_MEMORY_CACHE )
    ecbuild_add_library(
        TARGET atlas-orca
        TYPE SHARED
        INSTALL_HEADERS LISTED
        HEADER_DESTINATION "include/atlas-orca"
        PUBLIC_DEFINITIONS
            HAVE_ORCA_IN_MEMORY_CACHE
        SOURCES 
            Library.cc
            Library.h
            grid/Orca.cc
            grid/Orca.h
            grid/OrcaGrid.cc
            grid/OrcaGrid.h
            grid/OrcaGridBuilder.cc
            meshgenerator/OrcaMeshGenerator.cc
            meshgenerator/OrcaMeshGenerator.h
            util/AtlasIOReader.h
            util/ComputeCachedPath.h
            util/ComputeUid.h
            util/ComputeUid.cc
            util/DetectInvalidElements.cc
            util/DetectInvalidElements.h
            util/Download.cc
            util/Download.h
            util/Enums.h
            util/Flag.h
            util/OrcaData.cc
            util/OrcaData.h
            util/OrcaDataFile.h
            util/OrcaPeriodicity.cc
            util/OrcaPeriodicity.h
            util/PointIJ.h
            ${CMAKE_CURRENT_BINARY_DIR}/version.cc
            version.h
        PUBLIC_LIBS atlas
        PUBLIC_INCLUDES 
           $<INSTALL_INTERFACE:include/atlas-orca>
           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src> )
else()
    ecbuild_add_library(
        TARGET atlas-orca
        TYPE SHARED
        INSTALL_HEADERS LISTED
        HEADER_DESTINATION "include/atlas-orca"
        SOURCES 
            Library.cc
            Library.h
            grid/Orca.cc
            grid/Orca.h
            grid/OrcaGrid.cc
            grid/OrcaGrid.h
            grid/OrcaGridBuilder.cc
            meshgenerator/OrcaMeshGenerator.cc
            meshgenerator/OrcaMeshGenerator.h
            util/AtlasIOReader.h
            util/ComputeCachedPath.h
            util/ComputeUid.h
            util/ComputeUid.cc
            util/DetectInvalidElements.cc
            util/DetectInvalidElements.h
            util/Download.cc
            util/Download.h
            util/Enums.h
            util/Flag.h
            util/OrcaData.cc
            util/OrcaData.h
            util/OrcaDataFile.h
            util/OrcaPeriodicity.cc
            util/OrcaPeriodicity.h
            util/PointIJ.h
            ${CMAKE_CURRENT_BINARY_DIR}/version.cc
            version.h
        PUBLIC_LIBS atlas
        PUBLIC_INCLUDES 
           $<INSTALL_INTERFACE:include/atlas-orca>
           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src> )
endif()

if( atlas_VERSION VERSION_GREATER_EQUAL 0.35.0 )
    set( ATLAS_TEMPORARY_ELEMENTTYPES 0 )
elseif( NOT DEFINED ATLAS_TEMPORARY_ELEMENTTYPES )
    set( ATLAS_TEMPORARY_ELEMENTTYPES 1 )
endif()
set_source_files_properties( meshgenerator/OrcaMeshGenerator.cc PROPERTIES COMPILE_DEFINITIONS "ATLAS_TEMPORARY_ELEMENTTYPES=${ATLAS_TEMPORARY_ELEMENTTYPES}" )

if (eckit_HAVE_CURL)
    set_source_files_properties( util/Download.cc
                            PROPERTIES COMPILE_DEFINITIONS "ATLAS_ORCA_HAVE_ECKIT_URLHANDLE")
endif()
