# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
set(kernel-solver-list
  hypre.cpp
  superlu.cpp
  umfpack.cpp
  voxel_amavanka.cpp
  )

add_library(kernel-solver STATIC ${kernel-solver-list})
if (FEAT_COTIRE)
  cotire(kernel-solver)
endif (FEAT_COTIRE)
target_link_libraries(kernel-solver PRIVATE kernel-lafem kernel-util kernel-adjacency)

if (FEAT_HAVE_UMFPACK)
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/SuiteSparse_config/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/AMD/Include/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/UMFPACK/Include/")
endif (FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_HYPRE)
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_ls/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/utilities/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/seq_mv/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_mv/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/IJ_mv/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/multivector/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/krylov/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/distributed_matrix/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_block_mv/")
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/blas/")
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_SUPERLU_DIST)
  target_include_directories(kernel-solver PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUPERLU}/SRC/")
endif (FEAT_HAVE_SUPERLU_DIST)


if (FEAT_HAVE_CUDA)
  add_library(kernel-solver-cuda STATIC
    cusolver.cu
    ilu_precond.cu
    sor_precond.cu
    ssor_precond.cu
    voxel_amavanka.cu
  )
  target_link_libraries(kernel-solver-cuda PRIVATE ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cublas_LIBRARY} kernel-util-cuda)
  target_link_libraries(kernel-solver PRIVATE kernel-solver-cuda)
endif (FEAT_HAVE_CUDA)

# list of solver tests
set(test_list
  amavanka-test
  basic_solver-test
  cusolver-test
  hypre-test
  optimizer-test
  superlu-test
  umfpack-test
  vanka-test
  voxel_amavanka-test
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} feat test_system)

  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/solver/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    add_test(${test}_cuda_memcheck ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} cuda)
    set_property(TEST ${test}_cuda_memcheck PROPERTY LABELS "cuda memcheck")
    set_property(TEST ${test}_cuda_memcheck PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
endforeach (test)

# add all tests to lafem_tests
add_custom_target(solver_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests solver_tests)
