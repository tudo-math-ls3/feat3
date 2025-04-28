# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

#ADD_LIBRARY(kernel-space STATIC)
#TARGET_LINK_LIBRARIES (kernel-space kernel-util)

# list of space tests
set(test_list
  bogner_fox_schmit-test
  discontinuous-test
  element-regression-test
  lagrange1-test
  rannacher_turek-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET space_tests LINK_AGAINST feat-kernel-geometry-test-aux)
