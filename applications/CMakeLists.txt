# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_SUBDIRECTORY( meshopt_boundary-app )
ADD_SUBDIRECTORY( meshopt_r_adapt-app )
ADD_SUBDIRECTORY( meshopt_refinement-app )
ADD_SUBDIRECTORY( meshopt_screws-app )
ADD_SUBDIRECTORY( navier_stokes_screws-app )

# applications to build
SET (app_list
  foundation_psystem_test-app
  hierarch_transfer_test-app
  meshopt_hyperelasticity_resize-app
  navier_stokes_cp2d_q2q1
  navier_stokes_solver_factory
  poisson_dirichlet_2d
  poisson_neumann_2d
  poisson_solver_factory
  scarckin_poisson
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


if (FEAT_HAVE_MPI)
  ADD_TEST(sleep18 sleep 2)
  SET_PROPERTY(TEST sleep18 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(foundation_psystem_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_psystem_test-app
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_psystem_test-app ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST foundation_psystem_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_psystem_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(foundation_psystem_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_psystem_test-app
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/foundation_psystem_test-app)
  SET_PROPERTY(TEST foundation_psystem_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_psystem_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app ${MPIEXEC_POSTFLAGS})
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app ${MPIEXEC_POSTFLAGS})
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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/hierarch_transfer_test-app)
  SET_PROPERTY(TEST hierarch_transfer_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST hierarch_transfer_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### poisson_dirichlet_2d
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep1 sleep 2)
  SET_PROPERTY(TEST sleep1 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_2d_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep3 sleep 2)
  SET_PROPERTY(TEST sleep3 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_2d_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep4 sleep 2)
  SET_PROPERTY(TEST sleep4 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(poisson_dirichlet_2d_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test-iter 6 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_6 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_6 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(poisson_dirichlet_2d_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback)
  SET_PROPERTY(TEST poisson_dirichlet_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_dirichlet_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 21 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 31 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 66 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --test-iter 19 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback)
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

######################### scarckin_poisson
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep12 sleep 2)
  SET_PROPERTY(TEST sleep12 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(scarckin_poisson_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarckin_poisson
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarckin_poisson_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarckin_poisson_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep14 sleep 2)
  SET_PROPERTY(TEST sleep14 PROPERTY LABELS "mpi,sleep")

  #    if (FEAT_HAVE_CUDA)
  # ADD_TEST(scarckin_poisson_mpi_4_cuda ${CMAKE_CTEST_COMMAND}
  #   --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
  #   --build-generator ${CMAKE_GENERATOR}
  #   --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  #   --build-target scarckin_poisson
  #   --build-nocmake
  #   --build-noclean
  #   --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mem cuda --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  #SET_PROPERTY(TEST scarckin_poisson_mpi_4_cuda PROPERTY LABELS "mpi")
  # SET_PROPERTY(TEST scarckin_poisson_mpi_4_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  #endif(FEAT_HAVE_CUDA)

  ADD_TEST(scarckin_poisson_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarckin_poisson
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarckin_poisson_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarckin_poisson_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(sleep15 sleep 2)
  SET_PROPERTY(TEST sleep15 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(scarckin_poisson_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarckin_poisson
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST scarckin_poisson_mpi_6 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarckin_poisson_mpi_6 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

else (FEAT_HAVE_MPI)
  ADD_TEST(scarckin_poisson_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarckin_poisson
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback)
  SET_PROPERTY(TEST scarckin_poisson_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST scarckin_poisson_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(scarckin_poisson_serial_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target scarckin_poisson
      --build-nocmake
      --build-noclean
      --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/scarckin_poisson --level 5 0 --mem cuda --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback)
    SET_PROPERTY(TEST scarckin_poisson_serial_cuda PROPERTY LABELS "serial")
    SET_PROPERTY(TEST scarckin_poisson_serial_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_CUDA)
endif (FEAT_HAVE_MPI)

######################### navier_stokes_cp2d_q2q1

if (FEAT_HAVE_MPI)
  ADD_TEST(sleep30 sleep 2)
  SET_PROPERTY(TEST sleep30 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(navier_stokes_cp2d_q2q1_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_cp2d_q2q1
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_cp2d_q2q1 --setup nozzle --test-mode --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST navier_stokes_cp2d_q2q1_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST navier_stokes_cp2d_q2q1_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(navier_stokes_cp2d_q2q1_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_cp2d_q2q1
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/navier_stokes_cp2d_q2q1 --setup nozzle --test-mode --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/)
  SET_PROPERTY(TEST navier_stokes_cp2d_q2q1_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST navier_stokes_cp2d_q2q1_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
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
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 5 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini ${MPIEXEC_POSTFLAGS})
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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_solver_factory --test-iter 4 --level 5 0 --mesh ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml --parti-type manual 2level parmetis fallback --solver-ini ${FEAT_SOURCE_DIR}/data/ini/solver_example.ini)
  SET_PROPERTY(TEST poisson_solver_factory_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_solver_factory_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)


# add all tests to applications_tests
ADD_CUSTOM_TARGET(applications_tests DEPENDS ${app_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests applications_tests)
