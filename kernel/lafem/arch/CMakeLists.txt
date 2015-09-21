cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of arch source files
SET (kernel-lafem-arch-list
  axpy_generic.cpp
  axpy_generic_ell.cpp
  axpy_generic_banded.cpp
  component_invert_generic.cpp
  component_product_generic.cpp
  defect_generic.cpp
  defect_generic_ell.cpp
  defect_generic_banded.cpp
  diagonal_generic.cpp
  difference_generic.cpp
  dot_product_generic.cpp
  gather_prim_generic.cpp
  norm_generic.cpp
  product_matvec_generic.cpp
  product_matvec_generic_ell.cpp
  product_matvec_generic_banded.cpp
  scale_generic.cpp
  scale_row_col_generic.cpp
  slip_filter_generic.cpp
  sum_generic.cpp
  unit_filter_generic.cpp
  unit_filter_blocked_generic.cpp
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
if (FEAST_COTIRE)
  cotire(kernel-lafem-arch)
endif (FEAST_COTIRE)

if (FEAST_BACKENDS_MKL)
  TARGET_LINK_LIBRARIES(kernel-lafem-arch
     -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
  )
endif(FEAST_BACKENDS_MKL)

if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-arch-cuda
    axpy.cu
    component_invert.cu
    component_product.cu
    defect.cu
    difference.cu
    diagonal.cu
    dot_product.cu
    gather_prim.cu
    norm.cu
    product_matvec.cu
    scale.cu
    scale_row_col.cu
    slip_filter.cu
    sum.cu
    unit_filter.cu
    unit_filter_blocked.cu
    STATIC
  )
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-lafem-arch-cuda)
  target_link_libraries(kernel-lafem-arch-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-lafem-arch kernel-lafem-arch-cuda)
endif(FEAST_BACKENDS_CUDA)
