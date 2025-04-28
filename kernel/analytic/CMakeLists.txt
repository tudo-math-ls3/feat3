# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of solver tests
set(test_list
  auto_derive-test
  common_function-test
  lambda_function-test
  parsed_function-test
  distance_function-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET analytic_tests CUDA)
