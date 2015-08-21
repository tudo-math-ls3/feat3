cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
SET (kernel-solver-list
  umfpack.cpp
  )

if (FEAST_HAVE_UMFPACK)
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/SuiteSparse_config/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/AMD/Include/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/UMFPACK/Include/")
endif (FEAST_HAVE_UMFPACK)

ADD_LIBRARY(kernel-solver STATIC ${kernel-solver-list})
TARGET_LINK_LIBRARIES(kernel-solver kernel-lafem kernel-util kernel-adjacency)

# list of solver tests
SET (test_list
  basic_solver-test
  umfpack-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feast test_system)

  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/solver/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/solver/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_main PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_main PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)

  ADD_TEST(${test}_cuda ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${CUDA_MEMCHECK_EXE} ${FEAST_BINARY_DIR}/kernel/solver/${test} cuda cudadevicereset)
  SET_PROPERTY(TEST ${test}_cuda PROPERTY LABELS "cuda")
  if (FEAST_CUDAMEMCHECK AND FEAST_BACKENDS_CUDA)
    SET_PROPERTY(TEST ${test}_cuda PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAST_CUDAMEMCHECK AND FEAST_BACKENDS_CUDA)
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(solver_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests solver_tests)
