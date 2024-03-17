# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(navier_stokes_screws-app navier_stokes_screws-app.cpp)
target_link_libraries(navier_stokes_screws-app feat)

# With MPI, run the test in both debug and optimized builds
if (NOT FEAT_DEBUG_MODE)
  if (FEAT_HAVE_MPI)
    add_test(sleep1 sleep 2)
    set_property(TEST sleep1 PROPERTY LABELS "mpi,sleep")

    add_test(navier_stokes_screws_app_test_1_mpi_4 ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target navier_stokes_screws-app
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_screws-app/navier_stokes_screws-app --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes ${MPIEXEC_POSTFLAGS})
    set_property(TEST navier_stokes_screws_app_test_1_mpi_4 PROPERTY LABELS "mpi")
    set_property(TEST navier_stokes_screws_app_test_1_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  # Without MPI, run the test only in optimized builds
  else (FEAT_HAVE_MPI)
      add_test(navier_stokes_screws_app_test_1_serial ${CMAKE_CTEST_COMMAND}
        --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
        --build-generator ${CMAKE_GENERATOR}
        --build-makeprogram ${CMAKE_MAKE_PROGRAM}
        --build-target navier_stokes_screws-app
        --build-nocmake
        --build-noclean
        --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/navier_stokes_screws-app/navier_stokes_screws-app --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes)
      set_property(TEST navier_stokes_screws_app_test_1_serial PROPERTY LABELS "serial")
      set_property(TEST navier_stokes_screws_app_test_1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_MPI)
endif (NOT FEAT_DEBUG_MODE)
