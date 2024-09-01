# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of arch source files
set(kernel-lafem-arch-list
  empty_lafem_arch.cpp
  )

if (FEAT_EICKT)
  set(kernel-lafem-arch-list
    ${kernel-lafem-arch-list}
    axpy_generic-eickt.cpp
    apply_generic-eickt.cpp
    apply_generic_banded-eickt.cpp
    component_copy_generic-eickt.cpp
    component_invert_generic-eickt.cpp
    component_product_generic-eickt.cpp
    diagonal_generic-eickt.cpp
    dot_product_generic-eickt.cpp
    lumping_generic-eickt.cpp
    max_abs_index_generic-eickt.cpp
    max_index_generic-eickt.cpp
    min_abs_index_generic-eickt.cpp
    min_index_generic-eickt.cpp
    norm_generic-eickt.cpp
    product_matmat_generic-eickt.cpp
    row_norm_generic-eickt.cpp
    scale_generic-eickt.cpp
    scale_row_col_generic-eickt.cpp
    slip_filter_generic-eickt.cpp
    transpose_generic-eickt.cpp
    unit_filter_generic-eickt.cpp
    unit_filter_blocked_generic-eickt.cpp
    )
endif (FEAT_EICKT)

if (FEAT_HAVE_MKL)
  set(kernel-lafem-arch-list ${kernel-lafem-arch-list}
    axpy_mkl.cpp
    apply_mkl.cpp
    component_copy_mkl.cpp
    component_product_mkl.cpp
    dot_product_mkl.cpp
    max_abs_index_mkl.cpp
    min_abs_index_mkl.cpp
    norm_mkl.cpp
    product_matmat_mkl.cpp
    scale_mkl.cpp
    transpose_mkl.cpp
    unit_filter_mkl.cpp
  )
endif (FEAT_HAVE_MKL)


add_library(kernel-lafem-arch STATIC ${kernel-lafem-arch-list})
if (FEAT_COTIRE)
  cotire(kernel-lafem-arch)
endif (FEAT_COTIRE)

if (FEAT_HAVE_CUDA)
  add_library(kernel-lafem-arch-cuda STATIC
    axpy.cu
    apply.cu
    component_copy.cu
    component_invert.cu
    component_product.cu
    diagonal.cu
    dot_product.cu
    lumping.cu
    max_abs_index.cu
    min_abs_index.cu
    mirror.cu
    norm.cu
    product_matmat.cu
    row_norm.cu
    scale.cu
    scale_row_col.cu
    slip_filter.cu
    transpose.cu
    unit_filter.cu
    unit_filter_blocked.cu
  )

  target_link_libraries(kernel-lafem-arch-cuda PRIVATE kernel-util-cuda ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cublasLt_LIBRARY}) #with new cuda support, we want cmake to be able to link against stubs...
#TODO: Use Static variants?  target_link_libraries(kernel-lafem-arch-cuda kernel-util-cuda ${CUDA_cusparse_static_LIBRARY} ${CUDA_cublas_static_LIBRARY})
  target_link_libraries(kernel-lafem-arch PRIVATE kernel-lafem-arch-cuda)
endif (FEAT_HAVE_CUDA)
