
LIST( APPEND ATLAS_TEST_ENVIRONMENT
    ATLAS_ORCA_CACHING=1
    ATLAS_CACHE_PATH=${CMAKE_BINARY_DIR}/share )

include_directories( ${PROJECT_SOURCE_DIR}/src )

ecbuild_add_test( TARGET  atlas_test_orca_plugin
                  SOURCES test_orca_plugin.cc
                  LIBS    atlas
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

ecbuild_add_test( TARGET  atlas_test_orca_grid_specs
                  SOURCES test_orca_grid_specs.cc
                  LIBS    atlas
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

ecbuild_add_test( TARGET  atlas_test_orca_grid
                  SOURCES test_orca_grid.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4 )

ecbuild_add_test( TARGET  atlas_test_orca_mesh
                  SOURCES test_orca_mesh.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4 )

ecbuild_add_test( TARGET  atlas_test_orca_valid_elements
                  SOURCES test_orca_valid_elements.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4 )

ecbuild_add_test( TARGET  atlas_test_orca_polygon_locator
                  SOURCES test_orca_polygon_locator.cc
                  LIBS    atlas-orca
                  MPI     2
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4 AND eckit_HAVE_MPI)

ecbuild_add_test( TARGET  atlas_test_orca_mesh_boundaries_NOMPI
                  SOURCES test_orca_mesh_boundaries.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4)

ecbuild_add_test( TARGET  atlas_test_orca_mesh_boundaries_MPI2
                  SOURCES test_orca_mesh_boundaries.cc
                  LIBS    atlas-orca
                  MPI     2
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4  AND eckit_HAVE_MPI)

ecbuild_add_test( TARGET  atlas_test_orca_surrounding_rectangle_NOMPI
                  SOURCES test_orca_surrounding_rectangle.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4 )

ecbuild_add_test( TARGET  atlas_test_orca_surrounding_rectangle_MPI2
                  SOURCES test_orca_surrounding_rectangle.cc
                  LIBS    atlas-orca
                  MPI     2
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4  AND eckit_HAVE_MPI)

ecbuild_add_test( TARGET  atlas_test_orca_surrounding_rectangle_MPI4
                  SOURCES test_orca_surrounding_rectangle.cc
                  LIBS    atlas-orca
                  MPI     4
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4  AND eckit_HAVE_MPI)

ecbuild_add_test( TARGET  atlas_test_orca_local_grid_NOMPI
                  SOURCES test_orca_local_grid.cc
                  LIBS    atlas-orca
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4)

ecbuild_add_test( TARGET  atlas_test_orca_local_grid_MPI2
                  SOURCES test_orca_local_grid.cc
                  LIBS    atlas-orca
                  MPI     2
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  CONDITION eckit_HAVE_LZ4  AND eckit_HAVE_MPI)

