# CmakeLists for directory "benchmarks"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

set(benchmarks
  axpy-bench
  product_matvec-bench
)

ADD_CUSTOM_TARGET(benchmarks DEPENDS ${benchmarks})

FOREACH (bench ${benchmarks})
  ADD_EXECUTABLE(${bench} ${bench}.cpp)
  TARGET_LINK_LIBRARIES(${bench} feast)
ENDFOREACH(bench)
