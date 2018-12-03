# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_SUBDIRECTORY( meshopt_boundary-app )
ADD_SUBDIRECTORY( meshopt_r_adapt-app )
ADD_SUBDIRECTORY( meshopt_refinement-app )
ADD_SUBDIRECTORY( navier_stokes_ppnd )
ADD_SUBDIRECTORY( navier_stokes_screws-app )

# applications to build
SET (app_list
  hierarch_transfer_test-app
  meshopt_hyperelasticity_resize-app
  navier_stokes_ccnd_q2p1
  navier_stokes_cp2d_q2p1
  poisson_bench_mg
  poisson_dirichlet
  poisson_mixed
  poisson_neumann
  poisson_scarc
  poisson_solver_factory
  stokes_3field_bench
  stokes_dricav_2d
  stokes_poiseuille_2d
  stokes_vortex_2d
  stokes_solver_factory
  )

FOREACH (app ${app_list} )
  ADD_EXECUTABLE(${app} ${app}.cpp)
  TARGET_LINK_LIBRARIES(${app} feat)
  if(MPI_C_COMPILE_FLAGS)
    set_target_properties(${app} PROPERTIES
      COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_C_LINK_FLAGS)
    set_target_properties(${app} PROPERTIES
      LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()
ENDFOREACH(app)

#disable advanced variable tracking for debug informations on big applications (it would be disabled afterwards anyways
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"))
  set_source_files_properties(meshopt_hyperelasticity_resize-app.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
  set_source_files_properties(navier_stokes_ccnd_q2p1.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
  set_source_files_properties(navier_stokes_cp2d_q2pq.cpp APPEND_STRING PROPERTIES COMPILE_FLAGS " -fno-var-tracking-assignments")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"))

######################### hierarch_transfer_test-app

if (FEAT_HAVE_MPI)
  ADD_TEST(sleep317 sleep 2)
  SET_PROPERTY(TEST sleep317 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(hierarch_transfer_test-app_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target hierarch_transfer_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app 4 2:1 0 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST hierarch_transfer_test-app_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST hierarch_transfer_test-app_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep318 sleep 2)
  SET_PROPERTY(TEST sleep318 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(hierarch_transfer_test-app_mpi_16 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target hierarch_transfer_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app 6 4:1 2:1 0 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST hierarch_transfer_test-app_mpi_16 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST hierarch_transfer_test-app_mpi_16 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(hierarch_transfer_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target hierarch_transfer_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app 3 0)
  SET_PROPERTY(TEST hierarch_transfer_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST hierarch_transfer_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_dirichlet
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep1 sleep 2)
  SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep3 sleep 2)
  SET_PROPERTY(TEST sleep3 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep4 sleep 2)
  SET_PROPERTY(TEST sleep4 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet --test-iter 6 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_mpi_6 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_mpi_6 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(poisson_dirichlet_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST poisson_dirichlet_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_dirichlet_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_mixed
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep1 sleep 2)
  SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_mixed_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_mixed
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3
    ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_mixed --bc neumann --test-iter 8 --level 5 2 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_mixed_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_mixed_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(poisson_mixed_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_mixed
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE}
    ${FEAT_BINARY_DIR}/applications/poisson_mixed --bc neumann --test-iter 8 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST poisson_dirichlet_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_dirichlet_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### stokes_poiseuille_2d
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep5 sleep 2)
  SET_PROPERTY(TEST sleep5 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(stokes_poiseuille_2d_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 21 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep8 sleep 2)
  SET_PROPERTY(TEST sleep8 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(stokes_poiseuille_2d_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 31 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep19 sleep 2)
  SET_PROPERTY(TEST sleep19 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(stokes_poiseuille_2d_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 31 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_6 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_6 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(stokes_poiseuille_2d_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 19 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### navier_stokes_ccnd_q2p1

if (FEAT_HAVE_MPI)
  ADD_TEST(sleep30 sleep 2)
  SET_PROPERTY(TEST sleep30 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(navier_stokes_ccnd_q2p1_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_ccnd_q2p1
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_ccnd_q2p1 --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 3 0:1 0 --test-mode ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST navier_stokes_ccnd_q2p1_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST navier_stokes_ccnd_q2p1_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(navier_stokes_ccnd_q2p1_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_ccnd_q2p1
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/navier_stokes_ccnd_q2p1 --mesh ${FEAT_SOURCE_DIR}/data/meshes/flowbench_c2d_03_quad_64.xml --level 3 0 --test-mode)
  SET_PROPERTY(TEST navier_stokes_ccnd_q2p1_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST navier_stokes_ccnd_q2p1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### navier_stokes_cp2d_q2p1

if (FEAT_HAVE_MPI)
  ADD_TEST(sleep30 sleep 2)
  SET_PROPERTY(TEST sleep30 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(navier_stokes_cp2d_q2p1_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_cp2d_q2p1
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_cp2d_q2p1 --setup nozzle --test-mode --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST navier_stokes_cp2d_q2p1_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST navier_stokes_cp2d_q2p1_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(navier_stokes_cp2d_q2p1_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_cp2d_q2p1
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/navier_stokes_cp2d_q2p1 --setup nozzle --test-mode --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/)
  SET_PROPERTY(TEST navier_stokes_cp2d_q2p1_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST navier_stokes_cp2d_q2p1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_solver_factory
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep1 sleep 2)
  SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_solver_factory_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_solver_factory
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_solver_factory_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_solver_factory_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep31 sleep 2)
  SET_PROPERTY(TEST sleep31 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_solver_factory_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_solver_factory
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_solver_factory_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_solver_factory_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep32 sleep 2)
  SET_PROPERTY(TEST sleep32 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_solver_factory_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_solver_factory
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_solver_factory_mpi_6 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_solver_factory_mpi_6 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(poisson_solver_factory_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_solver_factory
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini)
  SET_PROPERTY(TEST poisson_solver_factory_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_solver_factory_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(poisson_solver_factory_serial_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target poisson_solver_factory
      --build-nocmake
      --build-noclean
      --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example_cuda.ini)
    SET_PROPERTY(TEST poisson_solver_factory_serial_cuda PROPERTY LABELS "serial")
    SET_PROPERTY(TEST poisson_solver_factory_serial_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_CUDA)
endif (FEAT_HAVE_MPI)


# add all tests to applications_tests
ADD_CUSTOM_TARGET(applications_tests DEPENDS ${app_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests applications_tests)
