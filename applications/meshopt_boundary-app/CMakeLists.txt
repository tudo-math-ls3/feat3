# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(meshopt_boundary-app meshopt_boundary-app.cpp)
target_link_libraries(meshopt_boundary-app feat)
add_dependencies(tests meshopt_boundary-app)

if (FEAT_HAVE_MPI)
  add_test(NAME meshopt_boundary_app_test_1_mpi_3
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:meshopt_boundary-app> --test 1 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(meshopt_boundary_app_test_1_mpi_3 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME meshopt_boundary_app_test_2_mpi_4
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:meshopt_boundary-app> --test 2 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(meshopt_boundary_app_test_2_mpi_4 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(NAME meshopt_boundary_app_test_1_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:meshopt_boundary-app> --test 1 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes")
  set_tests_properties(meshopt_boundary_app_test_1_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME meshopt_boundary_app_test_2_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:meshopt_boundary-app> --test 2 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes")
  set_tests_properties(meshopt_boundary_app_test_2_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
