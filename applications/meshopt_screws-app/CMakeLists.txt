# CmakeLists for directory "applications/meshopt_boundary-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(meshopt_screws-app meshopt_screws-app)
TARGET_LINK_LIBRARIES(meshopt_screws-app feat)

if (FEAT_HAVE_MPI)
  ADD_TEST(sleep1 sleep 2)
  SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(meshopt_screws_app_test_1_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_screws-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_screws-app/meshopt_screws-app --test 1 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST meshopt_screws_app_test_1_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST meshopt_screws_app_test_1_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep2 sleep 2)
  SET_PROPERTY(TEST sleep2 PROPERTY LABELS "mpi,sleep")

  if (NOT FEAT_DEBUG_MODE)
    ADD_TEST(meshopt_screws_app_test_2_mpi_4 ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target meshopt_screws-app
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/meshopt_screws-app/meshopt_screws-app --test 2 ${MPIEXEC_POSTFLAGS})
    SET_PROPERTY(TEST meshopt_screws_app_test_2_mpi_4 PROPERTY LABELS "mpi")
    SET_PROPERTY(TEST meshopt_screws_app_test_2_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (NOT FEAT_DEBUG_MODE)
else (FEAT_HAVE_MPI)
  ADD_TEST(meshopt_screws_app_test_1_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target meshopt_screws-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_screws-app/meshopt_screws-app --test 1)
  SET_PROPERTY(TEST meshopt_screws_app_test_1_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST meshopt_screws_app_test_1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  if (NOT FEAT_DEBUG_MODE)
    ADD_TEST(meshopt_screws_app_test_2_serial ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target meshopt_screws-app
      --build-nocmake
      --build-noclean
      --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/meshopt_screws-app/meshopt_screws-app --test 2)
    SET_PROPERTY(TEST meshopt_screws_app_test_2_serial PROPERTY LABELS "serial")
    SET_PROPERTY(TEST meshopt_screws_app_test_2_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (NOT FEAT_DEBUG_MODE)
endif (FEAT_HAVE_MPI)
