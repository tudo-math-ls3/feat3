cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# list of lafem source files
SET (kernel-lafem-list memory_pool.cpp)

if (FEAST_BACKENDS_MKL)
  SET (kernel-lafem-list ${kernel-lafem-list}
    axpy_mkl.cpp
    component_product_mkl.cpp
    defect_mkl.cpp
    difference_mkl.cpp
    dot_product_mkl.cpp
    norm_mkl.cpp
    product_matvec_mkl.cpp
    scale_mkl.cpp sum_mkl.cpp
  )
endif(FEAST_BACKENDS_MKL)

ADD_LIBRARY(kernel-lafem STATIC ${kernel-lafem-list})
if (FEAST_BACKENDS_MKL)
  TARGET_LINK_LIBRARIES(kernel-lafem
    libmkl_intel_ilp64.a
    libmkl_sequential.a
    libmkl_core.a
  )
endif(FEAST_BACKENDS_MKL)

if (FEAST_GMP)
TARGET_LINK_LIBRARIES(kernel-lafem
  libgmpxx.so
  libgmp.so
  )
endif (FEAST_GMP)


if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-cuda STATIC
    memory_pool.cu
    axpy.cu
    defect.cu
    difference.cu
    dot_product.cu
    component_product.cu
    norm.cu
    product_matvec.cu
    scale.cu
    sum.cu
  )
  ADD_DEPENDENCIES(kernel-lafem kernel-lafem-cuda)
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-lafem-cuda)
endif(FEAST_BACKENDS_CUDA)

# list of lafem tests
SET (test_list
  axpy-test
  bicgstab-test
  container_pool-test
  defect-test
  dense_matrix-test
  dense_vector-test
  difference-test
  dot_product-test
  component_product-test
  norm-test
  product-test
  product_matvec-test
  richardson-test
  scale-test
  sparse_matrix_conversion-test
  sparse_matrix_coo-test
  sparse_matrix_csr-test
  sparse_matrix_ell-test
  sum-test
  transposition-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feast test_system)

  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")

  ADD_TEST(${test}_cuda ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/lafem/${test}" cuda)
  SET_PROPERTY(TEST ${test}_cuda PROPERTY LABELS "cuda")
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(lafem_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests lafem_tests)
