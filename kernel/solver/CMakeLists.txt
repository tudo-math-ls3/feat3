# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
SET (kernel-solver-list
  hypre.cpp
  superlu.cpp
  umfpack.cpp
  )

ADD_LIBRARY(kernel-solver STATIC ${kernel-solver-list})
if (FEAT_COTIRE)
  cotire(kernel-solver)
endif (FEAT_COTIRE)
TARGET_LINK_LIBRARIES(kernel-solver kernel-lafem kernel-util kernel-adjacency)

if (FEAT_HAVE_UMFPACK)
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/SuiteSparse_config/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/AMD/Include/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/UMFPACK/Include/")
endif (FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_HYPRE)
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_ls/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/utilities/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/seq_mv/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_mv/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/IJ_mv/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/multivector/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/krylov/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/distributed_matrix/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_block_mv/")
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/blas/")
endif(FEAT_HAVE_HYPRE)

if (FEAT_HAVE_SUPERLU_DIST)
  SET_PROPERTY(TARGET kernel-solver APPEND PROPERTY INCLUDE_DIRECTORIES "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUPERLU}/SRC/")
endif (FEAT_HAVE_SUPERLU_DIST)


if (FEAT_HAVE_CUDA)
  CUDA_ADD_LIBRARY(kernel-solver-cuda
    cusolver.cu
    ilu_precond.cu
    sor_precond.cu
    ssor_precond.cu
    STATIC
  )
  target_link_libraries(kernel-solver-cuda ${CUDA_cusolver_LIBRARY})
  target_link_libraries(kernel-solver-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-solver kernel-solver-cuda)
endif(FEAT_HAVE_CUDA)

# list of solver tests
SET (test_list
  amavanka-test
  basic_solver-test
  cusolver-test
  hypre-test
  optimizer-test
  superlu-test
  umfpack-test
  vanka-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/solver/${test})
  SET_PROPERTY(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    ADD_TEST(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} generic)
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    ADD_TEST(${test}_cuda_memcheck ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} cuda)
    SET_PROPERTY(TEST ${test}_cuda_memcheck PROPERTY LABELS "cuda memcheck")
    SET_PROPERTY(TEST ${test}_cuda_memcheck PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    SET_PROPERTY(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    SET_PROPERTY(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(solver_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests solver_tests)
