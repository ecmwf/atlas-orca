# Compiler bugs force us to disable some tests

if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI"  AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 16.8 )
  set( atlas_fctest_parametrisation_DISABLED TRUE )
  ecbuild_warn( "PGI Fortran compiler version tested up to 16.7 has a compiler bug "
                "that leads to a segfaults for atlas_fctest_parametrisation. "
                "Disabling this test..." )
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  ## ATLAS-170 workaround internal compiler error

      ## Other workaround that makes ecbuild_remove_fortran_flags work when variables are empty
      string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_CAPS )
      if( NOT CMAKE_Fortran_FLAGS )
          set( CMAKE_Fortran_FLAGS " " )
      endif()
      if( NOT CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS} )
          set( CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS} " " )
      endif()

  ecbuild_remove_fortran_flags( "-g" )
endif()


if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA )
  set( ATLAS_TEST_ENVIRONMENT "ATLAS_RUN_NGPUS=1" )
endif()

ecbuild_add_test( TARGET  atlas_test_orca_plugin
                  SOURCES test_orca_plugin.cc
                  LIBS    atlas
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                      ATLAS_PLUGINS=atlas-orca
                      ATLAS_PLUGINS_SEARCH_PATHS=${CMAKE_BINARY_DIR} )

ecbuild_add_test( TARGET  atlas_test_orca
                  SOURCES test_orca.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

ecbuild_add_test( TARGET  atlas_test_orca_valid_elements
                  SOURCES test_orca_valid_elements.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
