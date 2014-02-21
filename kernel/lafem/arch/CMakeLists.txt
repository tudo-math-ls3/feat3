cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# list of arch source files
SET (kernel-lafem-arch-list
  axpy_generic.cpp
  component_product_generic.cpp
  defect_generic.cpp
  difference_generic.cpp
  dot_product_generic.cpp
  norm_generic.cpp
  product_matvec_generic.cpp
  scale_generic.cpp
  sum_generic.cpp
  )

if (FEAST_BACKENDS_MKL)
  SET (kernel-lafem-arch-list ${kernel-lafem-arch-list}
    axpy_mkl.cpp
    component_product_mkl.cpp
    defect_mkl.cpp
    difference_mkl.cpp
    dot_product_mkl.cpp
    norm_mkl.cpp
    product_matvec_mkl.cpp
    scale_mkl.cpp
    sum_mkl.cpp
  )
endif(FEAST_BACKENDS_MKL)

ADD_LIBRARY(kernel-lafem-arch STATIC ${kernel-lafem-arch-list})

if (FEAST_BACKENDS_MKL)
  TARGET_LINK_LIBRARIES(kernel-lafem-arch
    libmkl_intel_ilp64.a
    libmkl_sequential.a
    libmkl_core.a
  )
endif(FEAST_BACKENDS_MKL)

if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-arch-cuda
    axpy.cu
    component_product.cu
    defect.cu
    difference.cu
    dot_product.cu
    norm.cu
    product_matvec.cu
    scale.cu
    sum.cu
    STATIC
  )
  ADD_DEPENDENCIES(kernel-lafem-arch kernel-lafem-arch-cuda)
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-lafem-arch-cuda)
  target_link_libraries(kernel-lafem-arch-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-lafem-arch kernel-lafem-arch-cuda)
endif(FEAST_BACKENDS_CUDA)
