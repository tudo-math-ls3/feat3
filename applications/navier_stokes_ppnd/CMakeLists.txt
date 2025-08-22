# CmakeLists for directory "applications/navier_stokes_ppnd"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(navier_stokes_ppnd navier_stokes_ppnd.cpp)
target_link_libraries(navier_stokes_ppnd feat)

if (FEAT_HAVE_MPI)
  add_test(NAME navier_stokes_ppnd_test_mpi_16
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navier_stokes_ppnd> --level 1 --setup fb-c2d-03 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ --time-steps 1600 --tol-rel-a 0.99 --test-mode ${MPIEXEC_POSTFLAGS})
  set_tests_properties(navier_stokes_ppnd_test_mpi_16 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")

  add_test(NAME navier_stokes_ppnd_test_fullfixpoint_mpi_16
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:navier_stokes_ppnd> --level 1 --setup fb-c2d-03 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ --time-steps 1600 --tol-rel-a 0.99 --test-mode --no-nonlinear --fix-steps 400 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(navier_stokes_ppnd_test_fullfixpoint_mpi_16 PROPERTIES LABELS "mpi" FAILED_REGULAR_EXPRESSION "FAILED")
endif (FEAT_HAVE_MPI)

# add all tests to applications_tests
add_custom_target(navier_stokes_ppnd_tests DEPENDS navier_stokes_ppnd)

# build all tests through top lvl target tests
add_dependencies(tests navier_stokes_ppnd_tests)
