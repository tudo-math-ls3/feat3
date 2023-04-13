# CmakeLists for directory "benchmarks"

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

set(benchmarks
  axpy-bench
  block_product_matvec-bench
  meshperm_assembly-bench
  dot_product-bench
  mixedprec_multigrid-bench
  multiprec_hierarch-bench
  product_matmat-bench
  product_matcsrmat-bench
  product_matcsrmat-file-bench
  product_matvec-bench
  product_matvec_dense-bench
)

ADD_CUSTOM_TARGET(benchmarks DEPENDS ${benchmarks})

FOREACH (bench ${benchmarks})
  ADD_EXECUTABLE(${bench} ${bench}.cpp)
  TARGET_LINK_LIBRARIES(${bench} feat)
ENDFOREACH(bench)
