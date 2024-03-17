# CMakeLists for directory "applications/dfg95_bench"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add steady-state ccnd solver application for bench1
add_executable(dfg95-ccnd-steady ccnd_steady.cpp)
target_link_libraries(dfg95-ccnd-steady feat)
if (MPI_C_COMPILE_FLAGS)
  set_target_properties(dfg95-ccnd-steady PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif ()
if (MPI_C_LINK_FLAGS)
  set_target_properties(dfg95-ccnd-steady PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif ()

# add unsteady ccnd solver application for bench2/bench3
add_executable(dfg95-ccnd-unsteady ccnd_unsteady.cpp)
target_link_libraries(dfg95-ccnd-unsteady feat)
if (MPI_C_COMPILE_FLAGS)
  set_target_properties(dfg95-ccnd-unsteady PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif ()
if (MPI_C_LINK_FLAGS)
  set_target_properties(dfg95-ccnd-unsteady PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif ()

# add prototypic unsteady cp2d solver application for bench2/bench3
add_executable(dfg95-cp2d-unsteady cp2d_unsteady.cpp)
target_link_libraries(dfg95-cp2d-unsteady feat)
if (MPI_C_COMPILE_FLAGS)
  set_target_properties(dfg95-cp2d-unsteady PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif ()
if (MPI_C_LINK_FLAGS)
  set_target_properties(dfg95-cp2d-unsteady PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif ()

# quadmath versions of applications
if (FEAT_HAVE_QUADMATH)
  add_executable(dfg95-ccnd-steady-qp ccnd_steady.cpp)
  target_link_libraries(dfg95-ccnd-steady-qp feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(dfg95-ccnd-steady-qp PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_USE_QUADMATH")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(dfg95-ccnd-steady-qp PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_USE_QUADMATH")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(dfg95-ccnd-steady-qp PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()
  add_executable(dfg95-ccnd-unsteady-qp ccnd_unsteady.cpp)
  target_link_libraries(dfg95-ccnd-unsteady-qp feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(dfg95-ccnd-unsteady-qp PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_USE_QUADMATH")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(dfg95-ccnd-unsteady-qp PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_USE_QUADMATH")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(dfg95-ccnd-unsteady-qp PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()
endif (FEAT_HAVE_QUADMATH)

###############################################################################

if (FEAT_HAVE_MPI)
  add_test(dfg95-ccnd-steady-test-mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target dfg95-ccnd-steady
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/dfg95_bench/dfg95-ccnd-steady --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0:1 0 --bench 1 --test-mode ${MPIEXEC_POSTFLAGS})
  set_property(TEST dfg95-ccnd-steady-test-mpi PROPERTY LABELS "mpi")
  set_property(TEST dfg95-ccnd-steady-test-mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(dfg95-ccnd-steady-test-serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target dfg95-ccnd-steady
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/dfg95_bench/dfg95-ccnd-steady --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0 --bench 1 --test-mode)
  set_property(TEST dfg95-ccnd-steady-test-serial PROPERTY LABELS "serial")
  set_property(TEST dfg95-ccnd-steady-test-serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

###############################################################################

if (FEAT_HAVE_MPI)
  add_test(dfg95-ccnd-unsteady-test-mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target dfg95-ccnd-unsteady
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/dfg95_bench/dfg95-ccnd-unsteady --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0:1 0 --bench 2 --test-mode ${MPIEXEC_POSTFLAGS})
  set_property(TEST dfg95-ccnd-unsteady-test-mpi PROPERTY LABELS "mpi")
  set_property(TEST dfg95-ccnd-unsteady-test-mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(dfg95-ccnd-unsteady-test-serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target dfg95-ccnd-unsteady
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/dfg95_bench/dfg95-ccnd-unsteady --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 2 0 --bench 2 --test-mode)
  set_property(TEST dfg95-ccnd-unsteady-test-serial PROPERTY LABELS "serial")
  set_property(TEST dfg95-ccnd-unsteady-test-serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)
