# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

set(kernel-cubature-list
  empty_cubature.cpp
)

if (FEAT_EICKT)
  set(kernel-cubature-list
    ${kernel-cubature-list}
    dynamic_factory-eickt.cpp
  )
endif (FEAT_EICKT)

add_library(feat-kernel-cubature STATIC ${kernel-cubature-list})

target_link_libraries(feat-kernel-cubature PRIVATE feat-kernel-util)

# list of cubature tests
set(test_list
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} feat test_system)

  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/cubature/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/cubature/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

endforeach (test)

# add all tests to cubature_tests
add_custom_target(cubature_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests cubature_tests)
