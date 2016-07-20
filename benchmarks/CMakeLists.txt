# CmakeLists for directory "benchmarks"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

set(benchmarks
  axpy-bench
  block_product_matvec-bench
  product_matmat-bench
  product_matvec-bench
)

ADD_CUSTOM_TARGET(benchmarks DEPENDS ${benchmarks})

FOREACH (bench ${benchmarks})
  ADD_EXECUTABLE(${bench} ${bench}.cpp)
  TARGET_LINK_LIBRARIES(${bench} feat)
ENDFOREACH(bench)
