cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

cmake_policy(SET CMP0048 NEW)
if( POLICY CMP0074 )
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(ecbuild 3.1 REQUIRED)

project(atlas-orca VERSION 0.0.0 LANGUAGES CXX)

find_package(atlas 0.21.0 REQUIRED)

atlas_create_plugin( orca LIBRARY atlas_orca )

add_subdirectory("share")
add_subdirectory("src")


ecbuild_install_project(NAME ${PROJECT_NAME})

