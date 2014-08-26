# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

if (MPI_CXX_INCLUDE_PATH)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

#top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# applications to build
SET (app_list
  foundation_system_test-app
  foundation_comm_test-app
  foundation_synch_test-app
  foundation_global_ops_convex_2procs_test-app
  foundation_global_ops_convex_4procs_test-app
  foundation_global_ops_concave_3procs_test-app
  foundation_global_ops_concave_4procs_test-app
  foundation_global_ops_synch_scal_test-app
  foundation_partitioning_4procs_test-app
  scarc_test-app
  scarc_1proc_test-app
  scarc_gateway_test-app
  poisson_test-app
  poisson_ms2_test-app
  proto_stokes
  )

FOREACH (app ${app_list} )
  ADD_EXECUTABLE(${app} ${app}.cpp)
  TARGET_LINK_LIBRARIES(${app} feast ${MPI_CXX_LIBRARIES})
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(${app} PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif()
  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(${app} PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()
ENDFOREACH(app)


if (FEAST_MPI)
  ADD_TEST(foundation_comm_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_comm_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_comm_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_comm_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_comm_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_comm_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_comm_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/foundation_comm_test-app)
  SET_PROPERTY(TEST foundation_comm_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_comm_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_system_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_system_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_system_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_system_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_system_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_system_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_system_test-app
    --build-nocmake
    --build-noclean
    --test-systemand ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/foundation_system_test-app)
  SET_PROPERTY(TEST foundation_system_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_system_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(scarc_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/scarc_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarc_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarc_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(scarc_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/scarc_test-app)
  SET_PROPERTY(TEST scarc_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST scarc_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(scarc_1proc_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_1proc_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/scarc_1proc_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarc_1proc_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarc_1proc_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(scarc_1proc_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_1proc_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/scarc_1proc_test-app)
  SET_PROPERTY(TEST scarc_1proc_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST scarc_1proc_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_synch_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_synch_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_synch_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_synch_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_synch_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_synch_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_synch_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/foundation_synch_test-app)
  SET_PROPERTY(TEST foundation_synch_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_synch_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(poisson_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/poisson_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(poisson_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/poisson_test-app)
  SET_PROPERTY(TEST poisson_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_global_ops_convex_2procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_convex_2procs_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_convex_2procs_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_global_ops_convex_2procs_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_convex_2procs_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_global_ops_convex_2procs_test-app")
  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_global_ops_convex_4procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_convex_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_convex_4procs_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_global_ops_convex_4procs_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_convex_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_global_ops_convex_4procs_test-app")
  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_global_ops_concave_3procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_concave_3procs_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_concave_3procs_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_global_ops_concave_3procs_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_concave_3procs_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_global_ops_concave_3procs_test-app")
  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_partitioning_4procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_partitioning_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_partitioning_4procs_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_partitioning_4procs_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_partitioning_4procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_partitioning_4procs_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_partitioning_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_partitioning_4procs_test-app")
  SET_PROPERTY(TEST foundation_partitioning_4procs_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_partitioning_4procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)


if (FEAST_MPI)
  ADD_TEST(scarc_gateway_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_gateway_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/scarc_gateway_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarc_gateway_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarc_gateway_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(scarc_gateway_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_gateway_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/scarc_gateway_test-app")
  SET_PROPERTY(TEST scarc_gateway_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST scarc_gateway_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_global_ops_concave_4procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_concave_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_concave_4procs_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_global_ops_concave_4procs_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_concave_4procs_test-app
    --build-nocmake
    --build-noclean
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_global_ops_concave_4procs_test-app")
  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(poisson_ms2_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_ms2_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/poisson_ms2_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_ms2_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_ms2_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(poisson_ms2_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_ms2_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/poisson_ms2_test-app)
  SET_PROPERTY(TEST poisson_ms2_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_ms2_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

if (FEAST_MPI)
  ADD_TEST(foundation_global_ops_synch_scal_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_synch_scal_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_synch_scal_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAST_MPI)
  ADD_TEST(foundation_global_ops_synch_scal_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_global_ops_synch_scal_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/applications/foundation_global_ops_synch_scal_test-app)
  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAST_MPI)

# add all tests to applications_tests
ADD_CUSTOM_TARGET(applications_tests DEPENDS scarc_test-app foundation_system_test-app)

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests applications_tests)
