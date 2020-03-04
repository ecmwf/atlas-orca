cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

find_package(ecbuild 3.0 REQUIRED)

project(atlas-data VERSION 0.0.0 LANGUAGES CXX)

find_package(eckit REQUIRED)

add_subdirectory("share")
add_subdirectory("src")

ecbuild_install_project(NAME ${PROJECT_NAME})

