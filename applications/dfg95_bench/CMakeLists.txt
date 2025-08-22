# CMakeLists for directory "applications/dfg95_bench"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add steady-state ccnd solver application for bench1
add_executable(dfg95-ccnd-steady ccnd_steady.cpp)
target_link_libraries(dfg95-ccnd-steady feat)
add_dependencies(tests dfg95-ccnd-steady)

# add unsteady ccnd solver application for bench2/bench3
add_executable(dfg95-ccnd-unsteady ccnd_unsteady.cpp)
target_link_libraries(dfg95-ccnd-unsteady feat)
add_dependencies(tests dfg95-ccnd-unsteady)

# add prototypic unsteady cp2d solver application for bench2/bench3
add_executable(dfg95-cp2d-unsteady cp2d_unsteady.cpp)
target_link_libraries(dfg95-cp2d-unsteady feat)

# quadmath versions of applications
if (FEAT_HAVE_QUADMATH)
  add_executable(dfg95-ccnd-steady-qp ccnd_steady.cpp)
  target_link_libraries(dfg95-ccnd-steady-qp feat)
  set_target_properties(dfg95-ccnd-steady-qp PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_USE_QUADMATH")
  add_executable(dfg95-ccnd-unsteady-qp ccnd_unsteady.cpp)
  target_link_libraries(dfg95-ccnd-unsteady-qp feat)
  set_target_properties(dfg95-ccnd-unsteady-qp PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_USE_QUADMATH")
endif (FEAT_HAVE_QUADMATH)

###############################################################################

if (FEAT_HAVE_MPI)
  add_test(NAME dfg95-ccnd-steady-test-mpi
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:dfg95-ccnd-steady> --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0:1 0 --bench 1 --test-mode ${MPIEXEC_POSTFLAGS})
  set_tests_properties(dfg95-ccnd-steady-test-mpi PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(NAME dfg95-ccnd-steady-test-serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:dfg95-ccnd-steady> --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0 --bench 1 --test-mode)
  set_tests_properties(dfg95-ccnd-steady-test-serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

###############################################################################

if (FEAT_HAVE_MPI)
  add_test(NAME dfg95-ccnd-unsteady-test-mpi
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:dfg95-ccnd-unsteady> --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0:1 0 --bench 2 --test-mode ${MPIEXEC_POSTFLAGS})
  set_tests_properties(dfg95-ccnd-unsteady-test-mpi PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(NAME dfg95-ccnd-unsteady-test-serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:dfg95-ccnd-unsteady> --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0 --bench 2 --test-mode)
  set_tests_properties(dfg95-ccnd-unsteady-test-serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
