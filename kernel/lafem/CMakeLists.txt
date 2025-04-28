# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( arch )

# list of lafem source files
set(kernel-lafem-list
  empty_lafem.cpp
  )

if (FEAT_EICKT)
  set(kernel-lafem-list
    ${kernel-lafem-list}
    dense_vector-eickt.cpp
    sparse_vector-eickt.cpp
    sparse_matrix_csr-eickt.cpp
    sparse_matrix_cscr-eickt.cpp
    )
endif (FEAT_EICKT)

add_library(feat-kernel-lafem STATIC ${kernel-lafem-list})
if (FEAT_COTIRE)
  cotire(feat-kernel-lafem)
endif (FEAT_COTIRE)

target_link_libraries(feat-kernel-lafem
PUBLIC
  feat-dependency
PRIVATE
  feat-kernel-lafem-arch feat-kernel-util feat-kernel-adjacency)

if (FEAT_HAVE_CUDA)
  target_link_libraries(feat-kernel-lafem PRIVATE feat-kernel-util-cuda)
endif (FEAT_HAVE_CUDA)

# list of lafem tests
set(test_list
  dense_matrix-test
  dense_vector-test
  dense_vector_blocked-test
  matrix_mirror_buffer-test
  matrix_mult-test
  mean_filter-test
  mean_filter_blocked-test
  meta_filter-test
  meta_matrix-apply-test
  meta_matrix-io-test
  meta_mirror-test
  meta_to_scalar-test
  meta_vector-axpy-test
  meta_vector-comp_prod-test
  meta_vector-comp_invert-test
  meta_vector-dot-norm2-test
  meta_vector-io-test
  meta_vector-scale-test
  pointstar_factory-test
  slip_filter-test
  sparse_matrix_conversion-test
  sparse_matrix_banded-test
  sparse_matrix_csr-test
  sparse_matrix_cscr-test
  sparse_matrix_bcsr-test
  sparse_matrix_factory-test
  sparse_vector-test
  sparse_vector_blocked-test
  test_matrix_factory-test
  unit_filter-test
  unit_filter_blocked-test
  vector_mirror-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET lafem_tests CUDA)
