# CmakeLists for directory "applications/meshopt_boundary-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(meshopt_r_adapt-app meshopt_r_adapt-app)
TARGET_LINK_LIBRARIES(meshopt_r_adapt-app feat)

#if (FEAT_HAVE_MPI)
#  ADD_TEST(meshopt_r_adapt_app_mpi_4 ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target meshopt_r_adapt-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --testmode ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST meshopt_r_adapt_app_mpi_4 PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST meshopt_r_adapt_app_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
ADD_TEST(meshopt_r_adapt_app_serial ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target meshopt_r_adapt-app
  --build-nocmake
  --build-noclean
  --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --testmode )
SET_PROPERTY(TEST meshopt_r_adapt_app_serial PROPERTY LABELS "serial")
SET_PROPERTY(TEST meshopt_r_adapt_app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)
