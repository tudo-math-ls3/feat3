# CMakeLists for directory "area51/ccnd_fiber"

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_EXECUTABLE(ccnd-fiber-time-stepping-test ccnd_time_stepping_test.cpp)
TARGET_LINK_LIBRARIES(ccnd-fiber-time-stepping-test feat)

ADD_EXECUTABLE(ccnd-joined-fiber-time-stepping-test ccnd_joined_time_stepping_test.cpp)
TARGET_LINK_LIBRARIES(ccnd-joined-fiber-time-stepping-test feat)

#disable advanced variable tracking for debug informations on big applications (it would be disabled afterwards anyways
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(ccnd_time_stepping_test.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

ADD_EXECUTABLE(tensor_operations-test EXCLUDE_FROM_ALL tensor_operations-test.cpp)
TARGET_LINK_LIBRARIES(tensor_operations-test feat test_system)

ADD_TEST(tensor_operations-test_generic ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target tensor_operations-test
  --build-nocmake
  --build-noclean
  --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/tensor_operations-test generic)
SET_PROPERTY(TEST tensor_operations-test_generic PROPERTY LABELS "generic")
if (FEAT_VALGRIND)
  SET_PROPERTY(TEST tensor_operations-test_generic PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
  SET_PROPERTY(TEST tensor_operations-test_generic PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_VALGRIND)

if (FEAT_HAVE_MPI)
ADD_TEST(sleep305 sleep 2)
SET_PROPERTY(TEST sleep305 PROPERTY LABELS "mpi,sleep")
ADD_TEST(ccnd-fiber-time-stepping-test-mpi-lagrange1 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-mpi-lagrange1 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-mpi-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep307 sleep 2)
  SET_PROPERTY(TEST sleep307 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(ccnd-fiber-time-stepping-test-mpi-lagrange2 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-mpi-lagrange2 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-mpi-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep308 sleep 2)
SET_PROPERTY(TEST sleep308 PROPERTY LABELS "mpi,sleep")
ADD_TEST(ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-joined-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-joined-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --joined-sol ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep309 sleep 2)
  SET_PROPERTY(TEST sleep309 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-joined-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-joined-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 --joined-sol ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
ADD_TEST(ccnd-fiber-time-stepping-test-serial-lagrange1 ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target ccnd-fiber-time-stepping-test
  --build-nocmake
  --build-noclean
  --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3)
SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-serial-lagrange1 PROPERTY LABELS "serial")
SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-serial-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

ADD_TEST(ccnd-fiber-time-stepping-test-serial-lagrange2 ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target ccnd-fiber-time-stepping-test
  --build-nocmake
  --build-noclean
  --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2)
SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-serial-lagrange2 PROPERTY LABELS "serial")
SET_PROPERTY(TEST ccnd-fiber-time-stepping-test-serial-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
