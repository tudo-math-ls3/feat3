cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of arch source files
SET (kernel-lafem-arch-list
  empty_lafem_arch.cpp
  )

if (FEAT_EICKT)
  SET (kernel-lafem-arch-list
    ${kernel-lafem-arch-list}
    axpy_generic-eickt.cpp
    axpy_generic_ell-eickt.cpp
    axpy_generic_banded-eickt.cpp
    component_invert_generic-eickt.cpp
    component_product_generic-eickt.cpp
    diagonal_generic-eickt.cpp
    difference_generic-eickt.cpp
    dot_product_generic-eickt.cpp
    lumping_generic-eickt.cpp
    max_element_generic-eickt.cpp
    norm_generic-eickt.cpp
    product_matmat_generic-eickt.cpp
    scale_generic-eickt.cpp
    scale_row_col_generic-eickt.cpp
    scatter_axpy_prim_generic-eickt.cpp
    slip_filter_generic-eickt.cpp
    sum_generic-eickt.cpp
    unit_filter_generic-eickt.cpp
    unit_filter_blocked_generic-eickt.cpp
    )
endif (FEAT_EICKT)

if (FEAT_HAVE_MKL)
  SET (kernel-lafem-arch-list ${kernel-lafem-arch-list}
    axpy_mkl.cpp
    component_product_mkl.cpp
    difference_mkl.cpp
    dot_product_mkl.cpp
    max_element_mkl.cpp
    norm_mkl.cpp
    product_matmat_mkl.cpp
    scale_mkl.cpp
    sum_mkl.cpp
  )
endif(FEAT_HAVE_MKL)


ADD_LIBRARY(kernel-lafem-arch STATIC ${kernel-lafem-arch-list})
if (FEAT_COTIRE)
  cotire(kernel-lafem-arch)
endif (FEAT_COTIRE)

if (FEAT_HAVE_MKL)
  TARGET_LINK_LIBRARIES(kernel-lafem-arch
     -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
  )
endif(FEAT_HAVE_MKL)

if (FEAT_HAVE_CUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-arch-cuda
    axpy.cu
    component_invert.cu
    component_product.cu
    difference.cu
    diagonal.cu
    dot_product.cu
    lumping.cu
    max_element.cu
    norm.cu
    product_matmat.cu
    scale.cu
    scale_row_col.cu
    scatter_axpy_prim.cu
    slip_filter.cu
    sum.cu
    unit_filter.cu
    unit_filter_blocked.cu
    STATIC
  )
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-lafem-arch-cuda)
  target_link_libraries(kernel-lafem-arch-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-lafem-arch kernel-lafem-arch-cuda)
endif(FEAT_HAVE_CUDA)
