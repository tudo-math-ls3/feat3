# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_SUBDIRECTORY( meshopt_boundary-app )
ADD_SUBDIRECTORY( meshopt_refinement-app )
ADD_SUBDIRECTORY( meshopt_screws-app )

# applications to build
SET (app_list
  #  foundation_global_ops_convex_2procs_test-app
  #  foundation_global_ops_convex_4procs_test-app
  #  foundation_global_ops_concave_3procs_test-app
  #  foundation_global_ops_concave_4procs_test-app
  #  foundation_global_ops_synch_scal_test-app
  foundation_psystem_test-app
  meshopt_hyperelasticity_resize-app
  parti_poisson_dirichlet_2d
  poisson_dirichlet_2d
  poisson_neumann_2d
  stokes_dricav_2d
  stokes_poiseuille_2d
  stokes_vortex_2d
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

#if (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_convex_2procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_convex_2procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_convex_2procs_test-app ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_mpi PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_convex_2procs_test-app_serial ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_convex_2procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command "${FEAT_BINARY_DIR}/applications/foundation_global_ops_convex_2procs_test-app")
#  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_serial PROPERTY LABELS "serial")
#  SET_PROPERTY(TEST foundation_global_ops_convex_2procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)

#if (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_convex_4procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_convex_4procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_convex_4procs_test-app ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_mpi PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_convex_4procs_test-app_serial ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_convex_4procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command "${FEAT_BINARY_DIR}/applications/foundation_global_ops_convex_4procs_test-app")
#  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_serial PROPERTY LABELS "serial")
#  SET_PROPERTY(TEST foundation_global_ops_convex_4procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)

#if (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_concave_3procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_concave_3procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_concave_3procs_test-app ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_mpi PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_concave_3procs_test-app_serial ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_concave_3procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command "${FEAT_BINARY_DIR}/applications/foundation_global_ops_concave_3procs_test-app")
#  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_serial PROPERTY LABELS "serial")
#  SET_PROPERTY(TEST foundation_global_ops_concave_3procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)

#if (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_concave_4procs_test-app_mpi ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_concave_4procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_concave_4procs_test-app ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_mpi PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_concave_4procs_test-app_serial ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_concave_4procs_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command "${FEAT_BINARY_DIR}/applications/foundation_global_ops_concave_4procs_test-app")
#  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_serial PROPERTY LABELS "serial")
#  SET_PROPERTY(TEST foundation_global_ops_concave_4procs_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)

#if (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_synch_scal_test-app_mpi ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_synch_scal_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_synch_scal_test-app ${MPIEXEC_POSTFLAGS})
#  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_mpi PROPERTY LABELS "mpi")
#  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#else (FEAT_HAVE_MPI)
#  ADD_TEST(foundation_global_ops_synch_scal_test-app_serial ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target foundation_global_ops_synch_scal_test-app
#    --build-nocmake
#    --build-noclean
#    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/foundation_global_ops_synch_scal_test-app)
#  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_serial PROPERTY LABELS "serial")
#  SET_PROPERTY(TEST foundation_global_ops_synch_scal_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
#endif (FEAT_HAVE_MPI)


######################### poisson_dirichlet_2d
if (FEAT_HAVE_MPI)
  ADD_TEST(poisson_dirichlet_2d_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 5 --level 5 0 --part_min_elems 10 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(poisson_dirichlet_2d_mpi_4_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target poisson_dirichlet_2d
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 5 --level 5 0 --mem cuda --part_min_elems 10 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
    SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4_cuda PROPERTY LABELS "mpi")
    SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_4_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif(FEAT_HAVE_CUDA)

  ADD_TEST(poisson_dirichlet_2d_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 9 --level 5 0 --part_min_elems 100 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_dirichlet_2d_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(poisson_dirichlet_2d_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 9 --level 5 0 --part_min_elems 250 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 4 --level 5 0 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST poisson_dirichlet_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_dirichlet_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(poisson_dirichlet_2d_serial_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target poisson_dirichlet_2d
      --build-nocmake
      --build-noclean
      --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/poisson_dirichlet_2d --test_iter 4 --level 5 0 --mem cuda --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
    SET_PROPERTY(TEST poisson_dirichlet_2d_serial_cuda PROPERTY LABELS "serial")
    SET_PROPERTY(TEST poisson_dirichlet_2d_serial_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_CUDA)
endif (FEAT_HAVE_MPI)

######################### stokes_poiseuille_2d
if (FEAT_HAVE_MPI)
  ADD_TEST(stokes_poiseuille_2d_mpi_4 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --part_min_elems 10 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(stokes_poiseuille_2d_mpi_4_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target stokes_poiseuille_2d
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --mem cuda --part_min_elems 10 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
    SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4_cuda PROPERTY LABELS "mpi")
    SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_4_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif(FEAT_HAVE_CUDA)

  ADD_TEST(stokes_poiseuille_2d_mpi_3 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --part_min_elems 100 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_3 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST stokes_poiseuille_2d_mpi_3 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(stokes_poiseuille_2d_mpi_6 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target stokes_poiseuille_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 6 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --part_min_elems 250 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST stokes_poiseuille_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  if (FEAT_HAVE_CUDA)
    ADD_TEST(stokes_poiseuille_2d_serial_cuda ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target stokes_poiseuille_2d
      --build-nocmake
      --build-noclean
      --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/stokes_poiseuille_2d --level 5 0 --mem cuda --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
    SET_PROPERTY(TEST stokes_poiseuille_2d_serial_cuda PROPERTY LABELS "serial")
    SET_PROPERTY(TEST stokes_poiseuille_2d_serial_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_CUDA)
endif (FEAT_HAVE_MPI)

###### parti_poisson

if (FEAT_HAVE_MPI)
  ADD_TEST(parti_poisson_dirichlet_2d_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target parti_poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/parti_poisson_dirichlet_2d --level 5 3 --part_min_elems 64 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST parti_poisson_dirichlet_2d_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST parti_poisson_dirichlet_2d_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (FEAT_HAVE_MPI)
  ADD_TEST(parti_poisson_dirichlet_2d_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target parti_poisson_dirichlet_2d
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/applications/parti_poisson_dirichlet_2d --level 5 3 --part_min_elems 64 --meshfile ${FEAT_SOURCE_DIR}/data/meshes/unit-square-quad.xml)
  SET_PROPERTY(TEST parti_poisson_dirichlet_2d_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST parti_poisson_dirichlet_2d_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

# add all tests to applications_tests
ADD_CUSTOM_TARGET(applications_tests DEPENDS ${app_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests applications_tests)
