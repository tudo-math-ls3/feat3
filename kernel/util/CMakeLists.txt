cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of util source files
ADD_LIBRARY(kernel-util STATIC
  exception.cpp
  mesh_streamer.cpp
  property_map.cpp
  runtime.cpp
)


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

if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-util-cuda
    cuda_util.cu
    memory_pool.cu
    STATIC
  )
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-util-cuda)
  TARGET_LINK_LIBRARIES(kernel-util kernel-util-cuda)
  target_link_libraries(kernel-util-cuda ${CUDA_cusparse_LIBRARY})
endif(FEAST_BACKENDS_CUDA)

# list of util tests
SET (test_list
  assertion-test
  binary_stream-test
  linear_algebra-test
  math-test
  mesh_streamer-test
  meta_math-test
  property_map-test
  random-test
  simple_arg_parser-test
  string-test
  string_mapped-test
  tiny_algebra-test
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
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/util/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)
ENDFOREACH(test)

# add all tests to util_tests
ADD_CUSTOM_TARGET(util_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests util_tests)
