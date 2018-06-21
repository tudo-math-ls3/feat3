cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
SET (kernel-solver-list
  umfpack.cpp
  )

if (FEAT_EICKT)
  SET (kernel-solver-list
    ${kernel-solver-list}
    matrix_stock-eickt.cpp
    solver_factory-eickt1.cpp
    solver_factory-eickt2.cpp
    )
endif(FEAT_EICKT)

if (FEAT_HAVE_MKL)
  SET (kernel-solver-list ${kernel-solver-list}
    spai_precond_mkl.cpp
  )
endif(FEAT_HAVE_MKL)

ADD_LIBRARY(kernel-solver STATIC ${kernel-solver-list})
if (FEAT_COTIRE)
  cotire(kernel-solver)
endif (FEAT_COTIRE)
TARGET_LINK_LIBRARIES(kernel-solver kernel-lafem kernel-util kernel-adjacency)

if (FEAT_HAVE_CUDA)
  CUDA_ADD_LIBRARY(kernel-solver-cuda
    cusolver.cu
    ilu_precond.cu
    spai_precond.cu
    sor_precond.cu
    ssor_precond.cu
    STATIC
  )
  if (FEAT_HAVE_CUSOLVER)
    if (CMAKE_MAJOR_VERSION LESS 3)
      find_library(CUDA_cusolver_LIBRARY NAMES cusolver PATHS "${CUDA_TOOLKIT_ROOT_DIR}/lib64")
    endif (CMAKE_MAJOR_VERSION LESS 3)
    target_link_libraries(kernel-solver-cuda ${CUDA_cusolver_LIBRARY})
  endif (FEAT_HAVE_CUSOLVER)
  target_link_libraries(kernel-solver-cuda ${CUDA_cusparse_LIBRARY})
  TARGET_LINK_LIBRARIES(kernel-solver kernel-solver-cuda)
endif(FEAT_HAVE_CUDA)

# list of solver tests
SET (test_list
  amg-test
  basic_solver-test
  cusolver-test
  optimiser-test
  umfpack-test
  vanka-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_main PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_main PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  ADD_TEST(${test}_cuda ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/solver/${test} cuda)
  SET_PROPERTY(TEST ${test}_cuda PROPERTY LABELS "cuda")
  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    SET_PROPERTY(TEST ${test}_cuda PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(solver_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests solver_tests)
