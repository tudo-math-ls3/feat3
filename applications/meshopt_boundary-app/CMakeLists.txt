# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(meshopt_boundary-app meshopt_boundary-app.cpp)
target_link_libraries(meshopt_boundary-app feat)

if (FEAT_HAVE_MPI)
  add_test(sleep1 sleep 2)
  set_property(TEST sleep1 PROPERTY LABELS "mpi")

  add_test(meshopt_boundary_app_test_1_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --test 1 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes" ${MPIEXEC_POSTFLAGS})
  set_property(TEST meshopt_boundary_app_test_1_mpi_3 PROPERTY LABELS "mpi")
  set_property(TEST meshopt_boundary_app_test_1_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(sleep2 sleep 2)
  set_property(TEST sleep2 PROPERTY LABELS "mpi")

  add_test(meshopt_boundary_app_test_2_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --test 2 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes" ${MPIEXEC_POSTFLAGS})
  set_property(TEST meshopt_boundary_app_test_2_mpi_4 PROPERTY LABELS "mpi")
  set_property(TEST meshopt_boundary_app_test_2_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(meshopt_boundary_app_test_1_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --test 1 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes")
  set_property(TEST meshopt_boundary_app_test_1_serial PROPERTY LABELS "serial")
  set_property(TEST meshopt_boundary_app_test_1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  add_test(meshopt_boundary_app_test_2_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_boundary-app
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_boundary-app/meshopt_boundary-app --test 2 --mesh-path "${FEAT_SOURCE_DIR}/data/meshes")
  set_property(TEST meshopt_boundary_app_test_2_serial PROPERTY LABELS "serial")
  set_property(TEST meshopt_boundary_app_test_2_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
