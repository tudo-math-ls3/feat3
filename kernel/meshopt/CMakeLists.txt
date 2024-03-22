# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of meshopt tests
SET (kernel-meshopt-list
    empty_meshopt.cpp
    )

if (FEAT_EICKT)
  set (kernel-meshopt-list
    ${kernel-meshopt-list}
      dudv_functional-eickt.cpp
      hyperelasticity_functional-eickt.cpp
      rumpf_functional-eickt.cpp
      rumpf_trafo-eickt.cpp
    )
endif(FEAT_EICKT)

ADD_LIBRARY(kernel-meshopt STATIC ${kernel-meshopt-list})

TARGET_LINK_LIBRARIES (kernel-meshopt kernel-cubature kernel-geometry kernel-lafem kernel-util)

# list of meshopt tests
SET (test_list
    hyperelasticity_functional-test
    )

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/meshopt/${test})
  SET_PROPERTY(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    ADD_TEST(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/meshopt/${test} generic)
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)
ENDFOREACH(test)

# add all tests to meshopt_tests
ADD_CUSTOM_TARGET(meshopt_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests meshopt_tests)
