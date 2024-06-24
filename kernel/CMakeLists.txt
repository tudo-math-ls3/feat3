# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# include kernel subdirectories
add_subdirectory( adjacency )
add_subdirectory( analytic )
add_subdirectory( assembly )
add_subdirectory( cubature )
add_subdirectory( geometry )
add_subdirectory( global )
add_subdirectory( lafem )
add_subdirectory( meshopt )
add_subdirectory( space )
add_subdirectory( solver )
add_subdirectory( trafo )
add_subdirectory( util )
add_subdirectory( voxel_assembly )

# list of kernel source files
set(KERNEL_LIBRARY_FILES
  backend.cpp
  runtime.cpp
)

# add kernel-root library
add_library(kernel-root STATIC ${KERNEL_LIBRARY_FILES})

# kernel-root and kernel-util have cyclic dependency
target_link_libraries(kernel-root PRIVATE kernel-util)
target_link_libraries(kernel-util PRIVATE kernel-root)

# add main feat library
add_library(feat INTERFACE)

# list of libraries to include in feat lib
target_link_libraries(feat INTERFACE
  kernel-adjacency
# kernel-assembly
  kernel-cubature
  kernel-geometry
  kernel-lafem
  kernel-meshopt
# kernel-space
  kernel-solver
# kernel-trafo
  kernel-util
  kernel-voxel-assembly
  kernel-root
)

if (FEAT_HAVE_ALGLIB)
  target_link_libraries(feat INTERFACE thirdparty-alglib)
endif(FEAT_HAVE_ALGLIB)

if (FEAT_HAVE_CGAL)
  target_link_libraries(feat INTERFACE ${Boost_LIBRARIES})
endif(FEAT_HAVE_CGAL)

if (FEAT_HAVE_DEATH_HANDLER)
  target_link_libraries(feat INTERFACE thirdparty-deathhandler)
endif (FEAT_HAVE_DEATH_HANDLER)

if (FEAT_HAVE_FPARSER)
  target_link_libraries(feat INTERFACE thirdparty-fparser)
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_PARMETIS)
  target_link_libraries(feat INTERFACE thirdparty-parmetis)
endif (FEAT_HAVE_PARMETIS)

if (FEAT_HAVE_HYPRE)
  target_link_libraries(feat INTERFACE thirdparty-hypre)
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_SUPERLU_DIST)
  target_link_libraries(feat INTERFACE thirdparty-superlu)
endif(FEAT_HAVE_SUPERLU_DIST)

if (FEAT_HAVE_TRIANGLE)
  target_link_libraries(feat INTERFACE thirdparty-triangle)
endif (FEAT_HAVE_TRIANGLE)

if (FEAT_HAVE_TRILINOS)
  #LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
  target_link_libraries(feat INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif (FEAT_HAVE_TRILINOS)

if (FEAT_HAVE_UMFPACK)
  target_link_libraries(feat INTERFACE thirdparty-umfpack)
endif(FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_ZFP)
  target_link_libraries(feat INTERFACE thirdparty-zfp)
endif (FEAT_HAVE_ZFP)

if (FEAT_HAVE_ZLIB)
  target_link_libraries(feat INTERFACE thirdparty-zlib)
endif (FEAT_HAVE_ZLIB)

if (FEAT_HAVE_ZOLTAN)
  target_link_libraries(feat INTERFACE thirdparty-zoltan)
endif (FEAT_HAVE_ZOLTAN)

#libquadmath needs to be the last library linked against
if (FEAT_HAVE_QUADMATH)
  target_link_libraries(feat INTERFACE -lquadmath)
endif (FEAT_HAVE_QUADMATH)

# list of kernel tests
set(test_list
  shape-test
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
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test})
  SET_PROPERTY(TEST ${test}_all PROPERTY LABELS "all")

if (FEAT_VALGRIND)
  add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_VALGRIND)
endforeach (test)

# add all tests to kernel_tests
add_custom_target(kernel_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests kernel_tests)
