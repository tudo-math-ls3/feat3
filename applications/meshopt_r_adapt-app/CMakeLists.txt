# CmakeLists for directory "applications/meshopt_boundary-app"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(meshopt_r_adapt-app meshopt_r_adapt-app.cpp)
target_link_libraries(meshopt_r_adapt-app feat)

if (FEAT_HAVE_MPI)
  add_test(sleep1 sleep 2)
  set_property(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  add_test(meshopt_r_adapt_app_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_r_adapt-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --test 1 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app" ${MPIEXEC_POSTFLAGS})
  set_property(TEST meshopt_r_adapt_app_mpi_3 PROPERTY LABELS "mpi")
  set_property(TEST meshopt_r_adapt_app_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(sleep2 sleep 2)
  set_property(TEST sleep2 PROPERTY LABELS "mpi,sleep")

  add_test(meshopt_r_adapt_app_alignment_mpi_5 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_r_adapt-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --test 2 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app" ${MPIEXEC_POSTFLAGS})
  set_property(TEST meshopt_r_adapt_app_alignment_mpi_5 PROPERTY LABELS "mpi")
  set_property(TEST meshopt_r_adapt_app_alignment_mpi_5 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
add_test(meshopt_r_adapt_app_serial ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target meshopt_r_adapt-app
  --build-nocmake
  --build-noclean
  --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --test 1 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app")
set_property(TEST meshopt_r_adapt_app_serial PROPERTY LABELS "serial")
set_property(TEST meshopt_r_adapt_app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
add_test(meshopt_r_adapt_app_alignment_serial ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target meshopt_r_adapt-app
  --build-nocmake
  --build-noclean
  --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_r_adapt-app/meshopt_r_adapt-app --test 2 --mesh-path "${FEAT_SOURCE_DIR}/applications/meshopt_r_adapt-app")
set_property(TEST meshopt_r_adapt_app_alignment_serial PROPERTY LABELS "serial")
set_property(TEST meshopt_r_adapt_app_alignment_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
