cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

#list of lafem source files
ADD_LIBRARY(kernel-lafem STATIC memory_pool.cpp)

if (FEAST_BACKENDS_CUDA)
  INCLUDE(FindCUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-cuda STATIC memory_pool.cu axpy.cu)
  ADD_DEPENDENCIES(kernel-lafem kernel-lafem-cuda)
  set(CUDA_NVCC_FLAGS "-arch=sm_13" CACHE STRING "Semi-colon delimit multiple arguments." FORCE)
endif(FEAST_BACKENDS_CUDA)


#list of lafem tests
SET (test_list dense_vector-test dense_matrix-test sparse_matrix_coo-test sparse_matrix_csr-test sum-test
  element_product-test difference-test axpy-test dot_product-test norm-test product-test defect-test scale-test richardson-test bicgstab-test)

#create all tests
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
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_cpu ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" cpu)
  SET_PROPERTY(TEST ${test}_cpu PROPERTY LABELS "cpu")

  ADD_TEST(${test}_gpu ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" gpu)
  SET_PROPERTY(TEST ${test}_cpu PROPERTY LABELS "gpu")
ENDFOREACH(test)

#add all tests to lafem_tests
ADD_CUSTOM_TARGET(lafem_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests lafem_tests)
