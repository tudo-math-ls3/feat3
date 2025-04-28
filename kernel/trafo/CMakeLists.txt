# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

#ADD_LIBRARY(kernel-trafo STATIC)
#TARGET_LINK_LIBRARIES (kernel-trafo kernel-util)

# list of trafo tests
set(test_list
  standard_trafo-test
  inverse_mapping-test
  isoparam_trafo-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET trafo_tests)
