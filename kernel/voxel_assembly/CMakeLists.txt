# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( arch )

# list of lafem source files
set(kernel-voxel-assembly-list
  empty_voxel.cpp
  )

add_library(kernel-voxel-assembly STATIC ${kernel-voxel-assembly-list})
if (FEAT_COTIRE)
  cotire(kernel-voxel-assembly)
endif (FEAT_COTIRE)

target_link_libraries(kernel-voxel-assembly PRIVATE kernel-voxel-assembly-arch kernel-util)

if (FEAT_HAVE_CUDA)
  target_link_libraries(kernel-voxel-assembly PRIVATE kernel-util-cuda)
endif (FEAT_HAVE_CUDA)

# list of lafem tests
set(test_list
  voxel_assembly-test
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
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/voxel_assembly/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/voxel_assembly/${test} generic mkl)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    add_test(${test}_cuda_memcheck ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/voxel_assembly/${test} cuda)
    set_property(TEST ${test}_cuda_memcheck PROPERTY LABELS "cuda_memcheck")
    set_property(TEST ${test}_cuda_memcheck PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
endforeach (test)

# add all tests to voxel_assembly_tests
add_custom_target(voxel_assembly_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests voxel_assembly_tests)