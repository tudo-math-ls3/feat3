# CmakeLists for directory "control"

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

#list of test_system tests
set(test_list checkpoint-test)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET control_tests)
