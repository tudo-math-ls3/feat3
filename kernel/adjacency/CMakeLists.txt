# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of adjacency source files
add_library(feat-kernel-adjacency STATIC
  coloring.cpp
  cuthill_mckee.cpp
  graph.cpp
  permutation.cpp
)
if (FEAT_COTIRE)
  cotire(feat-kernel-adjacency)
endif (FEAT_COTIRE)

# list of adjacency tests
set(test_list
  coloring-test
  cuthill_mckee-test
  graph-test
  permutation-test
)

target_link_libraries(feat-kernel-adjacency PRIVATE feat-kernel-util feat-dependency)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET adjacency_tests)
