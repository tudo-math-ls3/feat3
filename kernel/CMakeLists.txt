cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# include kernel subdirectories
ADD_SUBDIRECTORY( adjacency )
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
ADD_LIBRARY(feat STATIC archs.cpp)

# list of libraries to include in feat lib
TARGET_LINK_LIBRARIES (feat
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
)

if (FEAT_HAVE_ALGLIB)
  TARGET_LINK_LIBRARIES (feat thirdparty-alglib)
endif(FEAT_HAVE_ALGLIB)

if (FEAT_HAVE_FPARSER)
  TARGET_LINK_LIBRARIES (feat thirdparty-fparser)
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_UMFPACK)
  TARGET_LINK_LIBRARIES (feat thirdparty-umfpack)
endif(FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_PARMETIS)
  TARGET_LINK_LIBRARIES (feat thirdparty-parmetis)
endif (FEAT_HAVE_PARMETIS)

if (FEAT_HAVE_ZLIB)
  TARGET_LINK_LIBRARIES(feat thirdparty-zlib)
endif (FEAT_HAVE_ZLIB)

if (FEAT_HAVE_TRIANGLE)
  TARGET_LINK_LIBRARIES(feat thirdparty-triangle)
endif (FEAT_HAVE_TRIANGLE)

if (FEAT_HAVE_HYPRE)
  TARGET_LINK_LIBRARIES(feat thirdparty-hypre)
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_CGAL)
  target_link_libraries(feat ${Boost_LIBRARIES})
endif(FEAT_HAVE_CGAL)

if (FEAT_HAVE_DEATH_HANDLER)
  TARGET_LINK_LIBRARIES(feat thirdparty-deathhandler)
endif (FEAT_HAVE_DEATH_HANDLER)


#libquadmath needs to be the last library linked against
if (FEAT_HAVE_QUADMATH)
  TARGET_LINK_LIBRARIES (feat -lquadmath)
endif (FEAT_HAVE_QUADMATH)

# list of kernel tests
SET (test_list
  shape-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)
  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
ENDFOREACH(test)

# add all tests to kernel_tests
ADD_CUSTOM_TARGET(kernel_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests kernel_tests)
