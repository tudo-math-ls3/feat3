# CmakeLists for directory "test_system"

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_library(test_system STATIC test_system.cpp)

## why public here, but tests are linked explicitly against feat
# target_link_libraries(test_system PUBLIC feat)
target_link_libraries(test_system PRIVATE feat)

#list of test_system tests
set(test_list test_system-test)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET test_system_tests)
