# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of global tests
set(test_list
  alg_dof_parti-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET global_tests)
