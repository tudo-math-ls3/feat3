# CmakeLists for directory "applications/meshopt_boundary-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(navier_stokes_screws-app navier_stokes_screws-app)
TARGET_LINK_LIBRARIES(navier_stokes_screws-app feat)

# With MPI, run the test in both debug and optimized builds
if (NOT FEAT_DEBUG_MODE)
  if (FEAT_HAVE_MPI)
    ADD_TEST(sleep1 sleep 2)
    SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

    ADD_TEST(navier_stokes_screws_app_test_1_mpi_4 ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target navier_stokes_screws-app
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_screws-app/navier_stokes_screws-app --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes ${MPIEXEC_POSTFLAGS})
    SET_PROPERTY(TEST navier_stokes_screws_app_test_1_mpi_4 PROPERTY LABELS "mpi")
    SET_PROPERTY(TEST navier_stokes_screws_app_test_1_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  # Without MPI, run the test only in optimized builds
  else (FEAT_HAVE_MPI)
      ADD_TEST(navier_stokes_screws_app_test_1_serial ${CMAKE_CTEST_COMMAND}
        --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
        --build-generator ${CMAKE_GENERATOR}
        --build-makeprogram ${CMAKE_MAKE_PROGRAM}
        --build-target navier_stokes_screws-app
        --build-nocmake
        --build-noclean
        --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/navier_stokes_screws-app/navier_stokes_screws-app --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes)
      SET_PROPERTY(TEST navier_stokes_screws_app_test_1_serial PROPERTY LABELS "serial")
      SET_PROPERTY(TEST navier_stokes_screws_app_test_1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_MPI)
endif (NOT FEAT_DEBUG_MODE)
