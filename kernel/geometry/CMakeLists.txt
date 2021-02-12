cmake_minimum_required (VERSION 3.9)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# include test_aux subdirectory
ADD_SUBDIRECTORY( test_aux )

# list of geometry source files and depenencies
SET(kernel-geometry-list
  cgal.cpp
)

if (FEAT_EICKT)
  SET (kernel-geometry-list
    ${kernel-geometry-list}
    conformal_mesh-eickt.cpp
    mesh_part-eickt.cpp
    mesh_node-eickt.cpp
    mesh_file_reader-eickt.cpp
  )
endif(FEAT_EICKT)

ADD_LIBRARY(kernel-geometry STATIC ${kernel-geometry-list})

TARGET_LINK_LIBRARIES (kernel-geometry kernel-util)

# list of geometry tests
SET (test_list
  boundary_factory-test
  cgal-test
  hit_test_factory-test
  index_calculator-test
  mesh_node-test-conf-quad
  mesh_part-test
  shape_convert-test
  standard_refinery-test-conf-quad
  standard_refinery-test-conf-hexa
  standard_refinery-test-conf-tria
  standard_refinery-test-conf-tetra
  standard_refinery-test-struct-quad
  standard_refinery-test-struct-hexa
  struct_index_mapping-test
  parsed_hit_test_factory-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  ADD_DEPENDENCIES(${test} kernel-geometry-test-aux)
  TARGET_LINK_LIBRARIES(${test} feat test_system kernel-geometry-test-aux)
  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/geometry/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/geometry/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_main PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_main PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

ENDFOREACH(test)


#disabling warnings with pragmas does not always work with gcc
#https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  SET_SOURCE_FILES_PROPERTIES(cgal.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-undef")
  SET_SOURCE_FILES_PROPERTIES(cgal-test.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-undef")

  #silence weird glibcxx warning
  #cgal/CGAL-4.14/include/CGAL/AABB_tree.h:1138:20:
  #warning: argument 1 value ‘18446744073709551615’ exceeds maximum object size 9223372036854775807 [-Walloc-size-larger-than=]
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0)
    SET_SOURCE_FILES_PROPERTIES(cgal.cpp APPEND_STRING PROPERTY COMPILE_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    SET_SOURCE_FILES_PROPERTIES(cgal-test.cpp APPEND_STRING PROPERTY COMPILE_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    SET_SOURCE_FILES_PROPERTIES(cgal.cpp APPEND_STRING PROPERTY LINK_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    SET_SOURCE_FILES_PROPERTIES(cgal-test.cpp APPEND_STRING PROPERTY LINK_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
  endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0)
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

# add all tests to geometry_tests
ADD_CUSTOM_TARGET(geometry_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests geometry_tests)
