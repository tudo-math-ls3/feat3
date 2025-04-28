# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

#ADD_LIBRARY(kernel-assembly STATIC)
#TARGET_LINK_LIBRARIES (kernel-assembly kernel-util)

# list of assembly tests
set(test_list
  bilinear_operator-test
  discrete_evaluator-test
  error_computer-test
  grid_transfer-test
  grid_transfer-mass-test
  grid_transfer-tip-test
  interpolator-test
  jump_stabil-test
  linear_functional-test
  mean_filter_assembler-test
  rew_projector-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET assembly_tests)
