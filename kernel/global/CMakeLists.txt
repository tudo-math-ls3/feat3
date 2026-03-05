# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of global tests
set(test_list
  alg_dof_parti-test
  remote_lambda-test
  global_inverse_mapping-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET global_tests)

if(FEAT_HAVE_MPI)
  add_test(NAME remote_lambda-test_4
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:remote_lambda-test> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(remote_lambda-test_4 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")

  add_test(NAME global_inverse_mapping-test_2
    COMMAND ${MPIEXEC} ${FEAT_MPI_MAPPING_FLAG} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:global_inverse_mapping-test> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(global_inverse_mapping-test_2 PROPERTIES LABELS "mpi" FAIL_REGULAR_EXPRESSION "FAILED")
endif()
