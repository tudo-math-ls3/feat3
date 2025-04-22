# CmakeLists for directory "benchmarks"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

set(benchmarks
  axpy-bench
  block_product_matvec-bench
  cuda_synch_bench
  meshperm_assembly-bench
  dot_product-bench
  io_bench
  mixedprec_multigrid-bench
  multiprec_hierarch-bench
  parperf1-bench
  parperf2-bench
  pcg-bench
  poisson_multigrid-bench
  product_matmat-bench
  product_matcsrmat-bench
  product_matcsrmat-file-bench
  product_matvec-bench
  product_matvec_dense-bench
  saddle_point_assembly-bench
  stokes_multigrid-bench
)

add_custom_target(benchmarks DEPENDS ${benchmarks})

foreach (bench ${benchmarks})
  add_executable(${bench} ${bench}.cpp)
  target_link_libraries(${bench} feat)
endforeach (bench)
