# CMakeLists for directory "applications"

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( ccnd )
add_subdirectory( ccnd_simple )
add_subdirectory( dfg95_bench )
add_subdirectory( gendie )
add_subdirectory( meshopt_boundary-app )
add_subdirectory( meshopt_r_adapt-app )
add_subdirectory( meshopt_refinement-app )
add_subdirectory( navier_stokes_ppnd )
add_subdirectory( navier_stokes_screws-app )
add_subdirectory( navsto_k_e_bfs )

# applications to build
set(app_list
  hierarch_transfer_test-app
  meshopt_hyperelasticity_resize-app
  poisson_bench_mg
  poisson_dirichlet
  poisson_mixed
  poisson_neumann
  poisson_simple_blocked
  poisson_simple_scalar
  poisson_voxel_mg
  poisson_scarc
  stokes_3field_bench
  stokes_dricav_2d
  stokes_poiseuille_2d
  stokes_simple
  stokes_vortex_2d
  stokes_voxel_ypipe
  write_stnd_matrices
  )

foreach (app ${app_list} )
  add_executable(${app} ${app}.cpp)
  target_link_libraries(${app} feat)
endforeach (app)

#disable advanced variable tracking for debug informations on big applications (it would be disabled afterwards anyways
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(meshopt_hyperelasticity_resize-app.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
  set_source_files_properties(navier_stokes_ccnd_q2p1.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
  set_source_files_properties(navier_stokes_cp2d_q2pq.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

######################### dfg95_bench1 applications

add_executable(dfg95-bench1-2d-dp dfg95_bench1.cpp)
add_executable(dfg95-bench1-2d-dp-sp dfg95_bench1.cpp)
add_executable(dfg95-bench1-3d-dp dfg95_bench1.cpp)
add_executable(dfg95-bench1-3d-dp-sp dfg95_bench1.cpp)

target_compile_definitions(dfg95-bench1-2d-dp    PUBLIC FEAT_DFG95_BENCH1_DIM=2 FEAT_DFG95_BENCH1_DOUBLE)
target_compile_definitions(dfg95-bench1-2d-dp-sp PUBLIC FEAT_DFG95_BENCH1_DIM=2 FEAT_DFG95_BENCH1_DOUBLE_SINGLE)
target_compile_definitions(dfg95-bench1-3d-dp    PUBLIC FEAT_DFG95_BENCH1_DIM=3 FEAT_DFG95_BENCH1_DOUBLE)
target_compile_definitions(dfg95-bench1-3d-dp-sp PUBLIC FEAT_DFG95_BENCH1_DIM=3 FEAT_DFG95_BENCH1_DOUBLE_SINGLE)

set(dfg95_bench1_list dfg95-bench1-2d-dp dfg95-bench1-2d-dp-sp dfg95-bench1-3d-dp dfg95-bench1-3d-dp-sp)

if (FEAT_HAVE_QUADMATH)
  add_executable(dfg95-bench1-2d-qp dfg95_bench1.cpp)
  add_executable(dfg95-bench1-2d-qp-dp dfg95_bench1.cpp)
  add_executable(dfg95-bench1-3d-qp dfg95_bench1.cpp)
  add_executable(dfg95-bench1-3d-qp-dp dfg95_bench1.cpp)

  target_compile_definitions(dfg95-bench1-2d-qp    PUBLIC FEAT_DFG95_BENCH1_DIM=2 FEAT_DFG95_BENCH1_QUAD)
  target_compile_definitions(dfg95-bench1-2d-qp-dp PUBLIC FEAT_DFG95_BENCH1_DIM=2 FEAT_DFG95_BENCH1_QUAD_DOUBLE)
  target_compile_definitions(dfg95-bench1-3d-qp    PUBLIC FEAT_DFG95_BENCH1_DIM=3 FEAT_DFG95_BENCH1_QUAD)
  target_compile_definitions(dfg95-bench1-3d-qp-dp PUBLIC FEAT_DFG95_BENCH1_DIM=3 FEAT_DFG95_BENCH1_QUAD_DOUBLE)

  set(dfg95_bench1_list ${dfg95_bench1_list} dfg95-bench1-2d-qp dfg95-bench1-2d-qp-dp dfg95-bench1-3d-qp dfg95-bench1-3d-qp-dp)
endif (FEAT_HAVE_QUADMATH)

foreach (app ${dfg95_bench1_list} )
  target_link_libraries(${app} feat)
endforeach (app)


######################### hierarch_transfer_test-app

if (FEAT_HAVE_MPI)
  add_test(NAME hierarch_transfer_test-app_mpi_4
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:hierarch_transfer_test-app> 4 2:1 0 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(hierarch_transfer_test-app_mpi_4 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME hierarch_transfer_test-app_mpi_16
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:hierarch_transfer_test-app> 6 4:1 2:1 0 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(hierarch_transfer_test-app_mpi_16 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  add_test(NAME hierarch_transfer_test-app_serial COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:hierarch_transfer_test-app> 3 0)
  set_tests_properties(hierarch_transfer_test-app_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_dirichlet
if (FEAT_HAVE_MPI)
  add_test(NAME poisson_dirichlet_mpi_4
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:poisson_dirichlet> --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(poisson_dirichlet_mpi_4 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
  set_property(TEST poisson_dirichlet_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME poisson_dirichlet_mpi_3
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:poisson_dirichlet> --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(poisson_dirichlet_mpi_3 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME poisson_dirichlet_mpi_6
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:poisson_dirichlet> --test-iter 6 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(poisson_dirichlet_mpi_6 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  add_test(NAME poisson_dirichlet_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:poisson_dirichlet> --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  set_tests_properties(poisson_dirichlet_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_mixed

if (FEAT_HAVE_MPI)
  add_test(NAME poisson_mixed_mpi_3
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:poisson_mixed> --bc neumann --test-iter 8 --level 5 2 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(poisson_mixed_mpi_3 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  add_test(NAME poisson_mixed_serial
    COMMAND ${VALGRIND_EXE} $<TARGET_FILE:poisson_mixed> --bc neumann --test-iter 8 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  set_tests_properties(poisson_mixed_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### stokes_poiseuille_2d

if (FEAT_HAVE_MPI)
  add_test(NAME stokes_poiseuille_2d_mpi_4
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:stokes_poiseuille_2d> --test-iter 21 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(stokes_poiseuille_2d_mpi_4 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME stokes_poiseuille_2d_mpi_3
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:stokes_poiseuille_2d> --test-iter 31 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(stokes_poiseuille_2d_mpi_3 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME stokes_poiseuille_2d_mpi_6
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:stokes_poiseuille_2d> --test-iter 31 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  set_tests_properties(stokes_poiseuille_2d_mpi_6 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  add_test(NAME stokes_poiseuille_2d_serial
    COMMAND ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} $<TARGET_FILE:stokes_poiseuille_2d> --test-iter 19 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  set_tests_properties(stokes_poiseuille_2d_serial PROPERTIES LABELS "serial" FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

# add all tests to applications_tests
add_custom_target(applications_tests DEPENDS ${app_list})

# build all tests through top lvl target tests
add_dependencies(tests applications_tests)
