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
  find_package(UMFPACK 6.3.5 REQUIRED)
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

if (FEAT_HAVE_TRILINOS)
  target_link_libraries(kernel-solver PUBLIC TRILINOS_TPL)
endif (FEAT_HAVE_TRILINOS)

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
  uzawa_bfbt-test
  vanka-test
  voxel_amavanka-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET solver_tests CUDA)
