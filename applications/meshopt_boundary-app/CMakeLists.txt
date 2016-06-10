# CmakeLists for directory "applications/meshopt_boundary-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(meshopt_boundary-app meshopt_boundary-app)
TARGET_LINK_LIBRARIES(meshopt_boundary-app feat)

if (FEAT_HAVE_MPI)
  ADD_TEST(meshopt_boundary_app_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --testmode ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST meshopt_boundary_app_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST meshopt_boundary_app_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(meshopt_boundary_app ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --testmode )
  SET_PROPERTY(TEST meshopt_boundary_app PROPERTY LABELS "serial")
  SET_PROPERTY(TEST meshopt_boundary_app PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
