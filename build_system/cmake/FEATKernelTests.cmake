# Helper function for adding FEAT3 kernel tests
#
# Parameters
# TESTS: One-value keyword argument, required, must contain the name of a list of test files
# TARGET: One-value keyword argument, required, must contain the name of a custom target
# CUDA: Option, if set CUDA versions of the given tests will be generated
# LINK_AGAINST: Many-value keyword argument, optional, contains additional targets tests will be linked against
#
# Description:
# Will generate the following tests for each file in the given TESTS list
# - default test (labeled with "all")
# - valgrind test, if FEAT_VALGRIND is true (labeled with "valgrind")
# - cuda test, if CUDA option is set and FEAT_HAVE_CUDA is true (labeled with "cuda")
# - cuda_memcheck test, if CUDA option is set and FEAT_HAVE_CUDA and FEAT_CUDAMEMCHECK are true (labeled with "cuda_memcheck")
#
# Usage:
# Assuming a set of tests
#
# set(test_list
#   foo.cpp
#   bar.cpp
#   baz.cpp)
#
# call this function either as
#
# add_feat_kernel_tests(TESTS test_list TARGET test_collection)
# or as
# add_feat_kernel_tests(TESTS test_list TARGET test_collection CUDA).
#
# foo.cpp, bar.cpp and baz.cpp will be added as tests and test_collection will
# be a common target that depends on all tests in the list.
# If the test have additional dependencies (most commonly
# feat-kernel-geometry-test-aux), you can link against these  using the
# LINK_AGAINST argument like this:
#
# add_feat_kernel_tests(
#   TEST test_list
#   TARGET test_collection
#   LINK_AGAINST feat-kernel-geometry-test-aux feat-kernel-solver)
#
# Tests are always linked against feat and test_system.
function(add_feat_kernel_tests)
  cmake_parse_arguments(
    PARSE_ARGV 0
    PARSED_ARGS
    "CUDA"
    "TESTS;TARGET"
    "LINK_AGAINST"
  )

  if(NOT PARSED_ARGS_TESTS)
    message(FATAL_ERROR "add_feat_kernel_tests: TESTS is a required keyword!")
  endif()

  if(NOT PARSED_ARGS_TARGET)
    message(FATAL_ERROR "add_feat_kernel_tests: TARGET is a required keyword!")
  endif()

  foreach (test ${${PARSED_ARGS_TESTS}} )
    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test} PRIVATE feat test_system)

    foreach (target ${PARSED_ARGS_LINK_AGAINST} )
      add_dependencies(${test} ${target})
      target_link_libraries(${test} PRIVATE ${target})
    endforeach()

    add_test(NAME ${test}_all COMMAND ${FEAT_CTEST_RUNNER} $<TARGET_FILE:${test}>)
    set_tests_properties(${test}_all PROPERTIES LABELS "all")

    if (FEAT_VALGRIND)
      add_test(NAME ${test}_valgrind COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:${test}> generic)
      set_tests_properties(${test}_valgrind PROPERTIES LABELS "valgrind")
      set_tests_properties(${test}_valgrind PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
      set_tests_properties(${test}_valgrind PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED")
    endif (FEAT_VALGRIND)

    if(PARSED_ARGS_CUDA)
      if(FEAT_HAVE_CUDA)
        add_test(NAME ${test}_cuda COMMAND ${FEAT_CTEST_RUNNER} $<TARGET_FILE:${test}> cuda)
        set_tests_properties(${test}_cuda PROPERTIES LABELS "cuda")
      endif()

      if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
        add_test(NAME ${test}_cuda_memcheck COMMAND ${FEAT_CTEST_RUNNER} ${CUDA_MEMCHECK_EXE} $<TARGET_FILE:${test}> cuda)
        set_tests_properties(${test}_cuda_memcheck PROPERTIES LABELS "cuda_memcheck")
        set_tests_properties(${test}_cuda_memcheck PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
        set_tests_properties(${test}_cuda_memcheck PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED")
        set_tests_properties(${test}_cuda_memcheck PROPERTIES FAIL_REGULAR_EXPRESSION "= Leaked")
      endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    endif(PARSED_ARGS_CUDA)

  endforeach (test)

  add_custom_target(${PARSED_ARGS_TARGET} DEPENDS ${${PARSED_ARGS_TESTS}})

  add_dependencies(tests ${PARSED_ARGS_TARGET})

  if(PARSED_ARGS_CUDA AND FEAT_HAVE_CUDA)
    add_dependencies(cuda_tests ${PARSED_ARGS_TARGET})
  endif()
endfunction()
