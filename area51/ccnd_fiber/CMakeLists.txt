# CMakeLists for directory "area51/ccnd_fiber"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(ccnd-fiber-time-stepping-test ccnd_time_stepping_test.cpp)
target_link_libraries(ccnd-fiber-time-stepping-test PRIVATE feat)
add_dependencies(tests ccnd-fiber-time-stepping-test)

add_executable(ccnd-joined-fiber-time-stepping-test ccnd_joined_time_stepping_test.cpp)
target_link_libraries(ccnd-joined-fiber-time-stepping-test PRIVATE feat)
add_dependencies(tests ccnd-joined-fiber-time-stepping-test)

#disable advanced variable tracking for debug informations on big applications (it would be disabled afterwards anyways
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(ccnd_time_stepping_test.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

add_executable(tensor_operations-test EXCLUDE_FROM_ALL tensor_operations-test.cpp)
target_link_libraries(tensor_operations-test PRIVATE feat test_system)
add_dependencies(tests tensor_operations-test)

add_test(NAME tensor_operations-test_all
  COMMAND ${FEAT_CTEST_RUNNER} $<TARGET_FILE:tensor_operations-test>)
set_property(TEST tensor_operations-test_all PROPERTY LABELS "all")

if (FEAT_VALGRIND)
  add_test(NAME tensor_operations-test_valgrind
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:tensor_operations-test> generic)
  set_property(TEST tensor_operations-test_valgrind PROPERTY LABELS "valgrind")
  set_property(TEST tensor_operations-test_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
  set_property(TEST tensor_operations-test_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_VALGRIND)

if (FEAT_HAVE_MPI)
  add_test(NAME ccnd-fiber-time-stepping-test-mpi-lagrange1
      COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ccnd-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ccnd-fiber-time-stepping-test-mpi-lagrange1 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

  add_test(NAME ccnd-fiber-time-stepping-test-mpi-lagrange2
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ccnd-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ccnd-fiber-time-stepping-test-mpi-lagrange2 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

  add_test(NAME ccnd-joined-fiber-time-stepping-test-mpi-lagrange1
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ccnd-joined-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --joined-sol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ccnd-joined-fiber-time-stepping-test-mpi-lagrange1 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

  add_test(NAME ccnd-joined-fiber-time-stepping-test-mpi-lagrange2
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:ccnd-joined-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 2 0:1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2 --joined-sol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(ccnd-joined-fiber-time-stepping-test-mpi-lagrange2 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  add_test(NAME ccnd-fiber-time-stepping-test-serial-lagrange1
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:ccnd-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3)
  set_tests_properties(ccnd-fiber-time-stepping-test-serial-lagrange1 PROPERTIES LABELS "serial" FAILED_REGULAR_EXPRESSION "FAILED")

  add_test(NAME ccnd-fiber-time-stepping-test-serial-lagrange2
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:ccnd-fiber-time-stepping-test> --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-cube-hexa.xml --level 1 0 --n-p 5.3 --n-s 3.9 --mu 0.1 --delta-t 0.0001 --testing-steps 3 --orient-space-2)
  set_tests_properties(ccnd-fiber-time-stepping-test-serial-lagrange2 PROPERTIES LABELS "serial" FAILED_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
