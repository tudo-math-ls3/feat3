# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(navier_stokes_screws-app navier_stokes_screws-app.cpp)
target_link_libraries(navier_stokes_screws-app feat)
add_dependencies(tests navier_stokes_screws-app)

# With MPI, run the test in both debug and optimized builds
if (NOT FEAT_DEBUG_MODE)
  if (FEAT_HAVE_MPI)
    add_test(NAME navier_stokes_screws_app_test_1_mpi_4
      COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navier_stokes_screws-app> --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes ${MPIEXEC_POSTFLAGS})
    set_tests_properties(navier_stokes_screws_app_test_1_mpi_4 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

  # Without MPI, run the test only in optimized builds
  else (FEAT_HAVE_MPI)
    add_test(NAME navier_stokes_screws_app_test_1_serial
      COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:navier_stokes_screws-app> --test 1 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes)
    set_tests_properties(navier_stokes_screws_app_test_1_serial PROPERTIES LABELS "serial" FAILED_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_MPI)
endif (NOT FEAT_DEBUG_MODE)
