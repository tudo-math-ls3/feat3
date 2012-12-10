cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

# list of geometry-test-aux source files and dependencies
ADD_LIBRARY(kernel-geometry-test-aux STATIC
  index_calculator_meshes.cpp
  standard_quad.cpp
  standard_hexa.cpp
  standard_tria.cpp
  standard_tetra.cpp
  tetris_quad.cpp
  tetris_hexa.cpp
  validate_structured_meshes.cpp
)

TARGET_LINK_LIBRARIES (kernel-geometry-test-aux kernel-util)
