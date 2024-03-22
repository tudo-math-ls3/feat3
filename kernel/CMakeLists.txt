# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# include kernel subdirectories
ADD_SUBDIRECTORY( adjacency )
ADD_SUBDIRECTORY( analytic )
ADD_SUBDIRECTORY( assembly )
ADD_SUBDIRECTORY( cubature )
ADD_SUBDIRECTORY( geometry )
ADD_SUBDIRECTORY( global )
ADD_SUBDIRECTORY( lafem )
ADD_SUBDIRECTORY( meshopt )
ADD_SUBDIRECTORY( space )
ADD_SUBDIRECTORY( solver )
ADD_SUBDIRECTORY( trafo )
ADD_SUBDIRECTORY( util )

# list of kernel source files
SET(KERNEL_LIBRARY_FILES
  backend.cpp
  runtime.cpp
)

# add kernel-root library
ADD_LIBRARY(kernel-root STATIC ${KERNEL_LIBRARY_FILES})

# kernel-root and kernel-util have cyclic dependency
TARGET_LINK_LIBRARIES(kernel-root kernel-util)
TARGET_LINK_LIBRARIES(kernel-util kernel-root)

# add main feat library
ADD_LIBRARY(feat INTERFACE)

# list of libraries to include in feat lib
TARGET_LINK_LIBRARIES (feat INTERFACE
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
  kernel-root
)

if (FEAT_HAVE_TRILINOS)
  #LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
  TARGET_LINK_LIBRARIES(feat INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif (FEAT_HAVE_TRILINOS)

if (FEAT_HAVE_ALGLIB)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-alglib)
endif(FEAT_HAVE_ALGLIB)

if (FEAT_HAVE_FPARSER)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-fparser)
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_UMFPACK)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-umfpack)
endif(FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_PARMETIS)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-parmetis)
endif (FEAT_HAVE_PARMETIS)

if (FEAT_HAVE_ZLIB)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-zlib)
endif (FEAT_HAVE_ZLIB)

if (FEAT_HAVE_TRIANGLE)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-triangle)
endif (FEAT_HAVE_TRIANGLE)

if (FEAT_HAVE_HYPRE)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-hypre)
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_SUPERLU_DIST)
  TARGET_LINK_LIBRARIES (feat INTERFACE thirdparty-superlu)
endif(FEAT_HAVE_SUPERLU_DIST)

if (FEAT_HAVE_CGAL)
  target_link_libraries(feat INTERFACE ${Boost_LIBRARIES})
endif(FEAT_HAVE_CGAL)

if (FEAT_HAVE_DEATH_HANDLER)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-deathhandler)
endif (FEAT_HAVE_DEATH_HANDLER)

if (FEAT_HAVE_ZFP)
  TARGET_LINK_LIBRARIES(feat INTERFACE thirdparty-zfp)
endif (FEAT_HAVE_ZFP)


#libquadmath needs to be the last library linked against
if (FEAT_HAVE_QUADMATH)
  TARGET_LINK_LIBRARIES(feat INTERFACE -lquadmath)
endif (FEAT_HAVE_QUADMATH)

# list of kernel tests
SET (test_list
  shape-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test})
  SET_PROPERTY(TEST ${test}_all PROPERTY LABELS "all")

if (FEAT_VALGRIND)
  ADD_TEST(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test} generic)
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif(FEAT_VALGRIND)
ENDFOREACH(test)

# add all tests to kernel_tests
ADD_CUSTOM_TARGET(kernel_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests kernel_tests)
