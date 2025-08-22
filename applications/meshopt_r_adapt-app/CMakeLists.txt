# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(meshopt_r_adapt-app meshopt_r_adapt-app.cpp)
target_link_libraries(meshopt_r_adapt-app feat)
add_dependencies(tests meshopt_r_adapt-app)

if (FEAT_HAVE_MPI)
  add_test(NAME meshopt_r_adapt_app_mpi_3
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:meshopt_r_adapt-app> --test 1 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(meshopt_r_adapt_app_mpi_3 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME meshopt_r_adapt_app_alignment_mpi_5
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:meshopt_r_adapt-app> --test 2 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(meshopt_r_adapt_app_alignment_mpi_5 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(NAME meshopt_r_adapt_app_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:meshopt_r_adapt-app> --test 1 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app")
  set_tests_properties(meshopt_r_adapt_app_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME meshopt_r_adapt_app_alignment_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:meshopt_r_adapt-app> --test 2 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app")
  set_tests_properties(meshopt_r_adapt_app_alignment_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
