# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
set(kernel-solver-list
  hypre.cpp
  mkl_dss.cpp
  superlu.cpp
  umfpack.cpp
  voxel_amavanka.cpp
  )

add_library(feat-kernel-solver STATIC ${kernel-solver-list})
if (FEAT_COTIRE)
  cotire(feat-kernel-solver)
endif (FEAT_COTIRE)
target_link_libraries(feat-kernel-solver
PUBLIC
  feat-dependency
PRIVATE
  feat-kernel-lafem feat-kernel-util feat-kernel-adjacency)

if (FEAT_HAVE_UMFPACK)
  find_package(UMFPACK 7.10.2 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE SuiteSparse::UMFPACK)
endif (FEAT_HAVE_UMFPACK)

if (FEAT_HAVE_HYPRE)
  find_package(HYPRE 2.33.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE HYPRE::HYPRE)
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_SUPERLU_DIST)
  find_package(superlu 9.1.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE SuperLU::SuperLU)
endif (FEAT_HAVE_SUPERLU_DIST)

if(FEAT_HAVE_ALGLIB)
  find_package(alglib 3.16.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PUBLIC Alglib::Alglib)
endif()


if (FEAT_HAVE_CUDA)
  add_library(feat-kernel-solver-cuda STATIC
    cudss.cu
    cusolver.cu
    ilu_precond.cu
    sor_precond.cu
    ssor_precond.cu
    voxel_amavanka.cu
  )
  if (FEAT_HAVE_CUDSS)
    target_link_libraries(feat-kernel-solver-cuda PRIVATE cuDSS::cuDSS)
  endif (FEAT_HAVE_CUDSS)
  target_link_libraries(feat-kernel-solver-cuda PRIVATE ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cublas_LIBRARY} feat-kernel-util-cuda)
  target_link_libraries(feat-kernel-solver PRIVATE feat-kernel-solver-cuda)
endif (FEAT_HAVE_CUDA)

# list of solver tests
set(test_list
  amavanka-test
  basic_solver-test
  cudss-test
  cusolver-test
  direct_stokes_solver-test
  hypre-test
  mkl_dss-test
  optimizer-test
  superlu-test
  umfpack-test
  vanka-test
  voxel_amavanka-test
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} PRIVATE feat test_system)

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
