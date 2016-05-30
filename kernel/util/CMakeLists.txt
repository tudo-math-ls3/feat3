cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of util source files
ADD_LIBRARY(kernel-util STATIC
  exception.cpp
  environment.cpp
  function_scheduler.cpp
  kahan_summation.cpp
  memory_pool.cpp
  memory_usage.cpp
  mpi_cout.cpp
  property_map.cpp
  runtime.cpp
  statistics.cpp
  xml_scanner.cpp
)

if (FEAT_COTIRE)
  cotire(kernel-util)
endif (FEAT_COTIRE)

#link against mpi if available
TARGET_LINK_LIBRARIES(kernel-util ${MPI_C_LIBRARIES})
if(MPI_C_COMPILE_FLAGS)
  set_target_properties(kernel-util PROPERTIES
    COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
if(MPI_C_LINK_FLAGS)
  set_target_properties(kernel-util PROPERTIES
    LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()

if (FEAT_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-util-cuda
    cuda_util.cu
    memory_pool.cu
    STATIC
  )
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-util-cuda)
  target_link_libraries(kernel-util-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-util kernel-util-cuda)
endif(FEAT_BACKENDS_CUDA)

# list of util tests
SET (test_list
  assertion-test
  binary_stream-test
  math-test
  memory_usage-test
  meta_math-test
  property_map-test
  random-test
  simple_arg_parser-test
  string-test
  string_mapped-test
  tiny_algebra-test
  xml_scanner-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)
  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/util/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)
ENDFOREACH(test)

# add all tests to util_tests
ADD_CUSTOM_TARGET(util_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests util_tests)
