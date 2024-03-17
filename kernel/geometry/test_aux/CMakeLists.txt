# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of geometry-test-aux source files and dependencies
add_library(kernel-geometry-test-aux STATIC
  index_calculator_meshes.cpp
  standard_quad.cpp
  standard_hexa.cpp
  standard_tria.cpp
  standard_tetra.cpp
  tetris_quad.cpp
  tetris_hexa.cpp
  validate_structured_meshes.cpp
)

target_link_libraries(kernel-geometry-test-aux PRIVATE kernel-util)
