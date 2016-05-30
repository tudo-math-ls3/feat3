cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# include test_aux subdirectory
ADD_SUBDIRECTORY( test_aux )

# list of geometry source files and depenencies
#ADD_LIBRARY(kernel-geometry STATIC)
#TARGET_LINK_LIBRARIES (kernel-geometry kernel-util)

# list of geometry tests
SET (test_list
  boundary_factory-test
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
ENDFOREACH(test)

# add all tests to geometry_tests
ADD_CUSTOM_TARGET(geometry_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests geometry_tests)
