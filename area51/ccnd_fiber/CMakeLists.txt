# CMakeLists for directory "area51/ccnd_fiber"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(ccnd-fiber-time-stepping-test ccnd_time_stepping_test.cpp)
target_link_libraries(ccnd-fiber-time-stepping-test feat)

add_executable(ccnd-joined-fiber-time-stepping-test ccnd_joined_time_stepping_test.cpp)
target_link_libraries(ccnd-joined-fiber-time-stepping-test feat)

#disable advanced variable tracking for debug informations on big applications (it would be disabled afterwards anyways
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(ccnd_time_stepping_test.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

add_executable(tensor_operations-test EXCLUDE_FROM_ALL tensor_operations-test.cpp)
target_link_libraries(tensor_operations-test feat test_system)

add_test(tensor_operations-test_all ${CMAKE_CTEST_COMMAND}
  --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-target tensor_operations-test
  --build-nocmake
  --build-noclean
  --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/tensor_operations-test)
set_property(TEST tensor_operations-test_all PROPERTY LABELS "all")

if (FEAT_VALGRIND)
  add_test(tensor_operations-test_valgrind ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target tensor_operations-test
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/tensor_operations-test generic)
  set_property(TEST tensor_operations-test_valgrind PROPERTY LABELS "valgrind")
  set_property(TEST tensor_operations-test_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
  set_property(TEST tensor_operations-test_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_VALGRIND)

if (FEAT_HAVE_MPI)
  add_test(sleep305 sleep 2)
  set_property(TEST sleep305 PROPERTY LABELS "mpi,sleep")
  add_test(ccnd-fiber-time-stepping-test-mpi-lagrange1 ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ccnd-fiber-time-stepping-test
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 ${MPIEXEC_POSTFLAGS})
  set_property(TEST ccnd-fiber-time-stepping-test-mpi-lagrange1 PROPERTY LABELS "mpi")
  set_property(TEST ccnd-fiber-time-stepping-test-mpi-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(sleep307 sleep 2)
  set_property(TEST sleep307 PROPERTY LABELS "mpi,sleep")

  add_test(ccnd-fiber-time-stepping-test-mpi-lagrange2 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 ${MPIEXEC_POSTFLAGS})
  set_property(TEST ccnd-fiber-time-stepping-test-mpi-lagrange2 PROPERTY LABELS "mpi")
  set_property(TEST ccnd-fiber-time-stepping-test-mpi-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(sleep308 sleep 2)
  set_property(TEST sleep308 PROPERTY LABELS "mpi,sleep")
  add_test(ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-joined-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-joined-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --joined-sol ${MPIEXEC_POSTFLAGS})
  set_property(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 PROPERTY LABELS "mpi")
  set_property(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(sleep309 sleep 2)
  set_property(TEST sleep309 PROPERTY LABELS "mpi,sleep")

  add_test(ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-joined-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-joined-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 --joined-sol ${MPIEXEC_POSTFLAGS})
  set_property(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 PROPERTY LABELS "mpi")
  set_property(TEST ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  add_test(ccnd-fiber-time-stepping-test-serial-lagrange1 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3)
  set_property(TEST ccnd-fiber-time-stepping-test-serial-lagrange1 PROPERTY LABELS "serial")
  set_property(TEST ccnd-fiber-time-stepping-test-serial-lagrange1 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(ccnd-fiber-time-stepping-test-serial-lagrange2 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ccnd-fiber-time-stepping-test
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/area51/ccnd_fiber/ccnd-fiber-time-stepping-test --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2)
  set_property(TEST ccnd-fiber-time-stepping-test-serial-lagrange2 PROPERTY LABELS "serial")
  set_property(TEST ccnd-fiber-time-stepping-test-serial-lagrange2 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
