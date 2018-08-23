# CmakeLists for directory "applications/navier_stokes_ppnd"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(navier_stokes_ppnd navier_stokes_ppnd)
TARGET_LINK_LIBRARIES(navier_stokes_ppnd feat)
  if(MPI_C_COMPILE_FLAGS)
    set_target_properties(navier_stokes_ppnd PROPERTIES
      COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif()
  if(MPI_C_LINK_FLAGS)
    set_target_properties(navier_stokes_ppnd PROPERTIES
      LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif()

  ADD_TEST(navier_stokes_ppnd_test_mpi_16 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_ppnd
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_ppnd/navier_stokes_ppnd --level 1 --setup fb-c2d-03 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ --time-steps 1600 --tol-rel-a 0.99 --test-mode ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST navier_stokes_ppnd_test_mpi_16 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST navier_stokes_ppnd_test_mpi_16 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

  ADD_TEST(navier_stokes_ppnd_test_fullfixpoint_mpi_16 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target navier_stokes_ppnd
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16 ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/applications/navier_stokes_ppnd/navier_stokes_ppnd --level 1 --setup fb-c2d-03 --mesh-path ${FEAT_SOURCE_DIR}/data/meshes/ --time-steps 1600 --tol-rel-a 0.99 --test-mode --no-nonlinear --fix-steps 400 ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST navier_stokes_ppnd_test_fullfixpoint_mpi_16 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST navier_stokes_ppnd_test_fullfixpoint_mpi_16 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")

# add all tests to applications_tests
ADD_CUSTOM_TARGET(navier_stokes_ppnd_tests DEPENDS navier_stokes_ppnd)

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests navier_stokes_ppnd_tests)
