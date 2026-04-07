# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of arch source files
set(kernel-lafem-arch-list
  empty_lafem_arch.cpp
  )

if (FEAT_EICKT)
  set(kernel-lafem-arch-list ${kernel-lafem-arch-list}
    axpy_dense_generic-eickt.cpp
    component_invert_dense_generic-eickt.cpp
    component_product_dense_generic-eickt.cpp
    copy_stride_generic-eickt.cpp
    diagonal_bcsr_generic-eickt.cpp
    diagonal_csr_generic-eickt.cpp
    dot_product_dense_generic-eickt.cpp
    lumping_bcsr_generic-eickt.cpp
    lumping_csr_generic-eickt.cpp
    matmatmult_dense_dense_generic-eickt.cpp
    matmatmult_sparse_dense_generic-eickt.cpp
    matvecmult_banded_dense_generic-eickt.cpp
    matvecmult_bcsr_dense_generic-eickt.cpp
    matvecmult_cscr_dense_generic-eickt.cpp
    matvecmult_csr_block_generic-eickt.cpp
    matvecmult_csr_dense_generic-eickt.cpp
    matvecmult_dense_dense_generic-eickt.cpp
    max_rel_diff_dense_generic-eickt.cpp
    min_max_index_dense_generic-eickt.cpp
    min_max_value_dense_generic-eickt.cpp
    norm2_dense_generic-eickt.cpp
    row_norm2_bcsr_generic-eickt.cpp
    row_norm2_csr_generic-eickt.cpp
    scale_col_bcsr_generic-eickt.cpp
    scale_col_csr_generic-eickt.cpp
    scale_dense_generic-eickt.cpp
    scale_row_bcsr_generic-eickt.cpp
    scale_row_csr_generic-eickt.cpp
    slip_filter_block_generic-eickt.cpp
    transpose_dense_generic-eickt.cpp
    triple_dot_product_dense_generic-eickt.cpp
    unit_filter_block_bcsr_generic-eickt.cpp
    unit_filter_block_vec_generic-eickt.cpp
    unit_filter_block_weak_bcsr_generic-eickt.cpp
    unit_filter_dense_csr_generic-eickt.cpp
    unit_filter_dense_vec_generic-eickt.cpp
    unit_filter_dense_weak_csr_generic-eickt.cpp
    )
endif (FEAT_EICKT)

if (FEAT_HAVE_MKL)
  set(kernel-lafem-arch-list ${kernel-lafem-arch-list}
    axpy_dense_mkl.cpp
    component_product_dense_mkl.cpp
    copy_stride_mkl.cpp
    dot_product_dense_mkl.cpp
    matmatmult_dense_dense_mkl.cpp
    matvecmult_bcsr_dense_mkl.cpp
    matvecmult_csr_dense_mkl.cpp
    matvecmult_dense_dense_mkl.cpp
    min_max_index_dense_mkl.cpp
    min_max_value_dense_mkl.cpp
    norm2_dense_mkl.cpp
    scale_dense_mkl.cpp
    transpose_dense_mkl.cpp
  )
endif (FEAT_HAVE_MKL)


add_library(feat-kernel-lafem-arch STATIC ${kernel-lafem-arch-list})
if (FEAT_COTIRE)
  cotire(feat-kernel-lafem-arch)
endif (FEAT_COTIRE)

target_link_libraries(feat-kernel-lafem-arch PUBLIC feat-dependency)

if (FEAT_HAVE_CUDA)
  add_library(feat-kernel-lafem-arch-cuda STATIC
    axpy_dense.cu
    component_invert_dense.cu
    component_product_dense.cu
    copy_stride.cu
    diagonal_bcsr.cu
    diagonal_csr.cu
    dot_product_dense.cu
    lumping_bcsr.cu
    lumping_csr.cu
    matmatmult_dense_dense.cu
    matmatmult_sparse_dense.cu
    matvecmult_banded_dense.cu
    matvecmult_bcsr_dense.cu
    matvecmult_csr_block.cu
    matvecmult_csr_dense.cu
    matvecmult_dense_dense.cu
    min_max_index_dense.cu
    min_max_value_dense.cu
    mirror_block_gather.cu
    mirror_block_scatter.cu
    mirror_dense_gather.cu
    mirror_dense_scatter.cu
    norm2_dense.cu
    row_norm2_bcsr.cu
    row_norm2_csr.cu
    scale_col_bcsr.cu
    scale_col_csr.cu
    scale_dense.cu
    scale_row_bcsr.cu
    scale_row_csr.cu
    slip_filter_block.cu
    transpose_dense.cu
    triple_dot_product_dense.cu
    unit_filter_block_bcsr.cu
    unit_filter_block_vec.cu
    unit_filter_block_weak_bcsr.cu
    unit_filter_dense_csr.cu
    unit_filter_dense_vec.cu
    unit_filter_dense_weak_csr.cu
  )

  target_link_libraries(feat-kernel-lafem-arch-cuda
  PUBLIC
    feat-dependency
  PRIVATE
    feat-kernel-util-cuda ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cublasLt_LIBRARY}) #with new cuda support, we want cmake to be able to link against stubs...
#TODO: Use Static variants?  target_link_libraries(kernel-lafem-arch-cuda kernel-util-cuda ${CUDA_cusparse_static_LIBRARY} ${CUDA_cublas_static_LIBRARY})

  if(FEAT_HAVE_CUDSS)
    find_package(cuDSS)
    target_link_libraries(feat-kernel-lafem-arch-cuda PRIVATE cuDSS::cuDSS)
  endif()

  target_link_libraries(feat-kernel-lafem-arch PRIVATE feat-kernel-lafem-arch-cuda)
endif (FEAT_HAVE_CUDA)
