# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of lafem source files
set(kernel-solver-list
  cudss.cpp
  direct_sparse_solver_cudss.cpp
  direct_sparse_solver_mkl.cpp
  direct_sparse_solver_mumps.cpp
  direct_sparse_solver_superlu.cpp
  direct_sparse_solver_umfpack.cpp
  frosch.cpp
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

if (FEAT_HAVE_CUDA)
  add_library(feat-kernel-solver-cuda STATIC
    ilu_precond.cu
    sor_precond.cu
    ssor_precond.cu
    voxel_amavanka.cu
  )

  target_link_libraries(feat-kernel-solver-cuda PRIVATE ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY} ${CUDA_cublas_LIBRARY} feat-kernel-util-cuda)
  target_link_libraries(feat-kernel-solver PRIVATE feat-kernel-solver-cuda)
endif (FEAT_HAVE_CUDA)

if(FEAT_HAVE_ALGLIB)
  find_package(Alglib 3.16.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PUBLIC Alglib::Alglib)
endif()

if(FEAT_HAVE_CUDSS)
  target_link_libraries(feat-kernel-solver PRIVATE cuDSS::cuDSS)
  target_include_directories(feat-kernel-solver PRIVATE "${CUDAToolkit_INCLUDE_DIRS}")
endif()

if (FEAT_HAVE_HYPRE)
  find_package(HYPRE 2.33.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE HYPRE::HYPRE)
endif (FEAT_HAVE_HYPRE)

if (FEAT_HAVE_MUMPS)
  target_link_libraries(feat-kernel-solver PRIVATE mumps::mumps)
endif (FEAT_HAVE_MUMPS)

if (FEAT_HAVE_SUPERLU_DIST)
  find_package(superlu 9.1.0 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE SuperLU::SuperLU)
endif (FEAT_HAVE_SUPERLU_DIST)

if (FEAT_HAVE_TRILINOS)
  find_package(Trilinos
      REQUIRED
      PATHS "${FEAT_TRILINOS_PATH}/"
      COMPONENTS
      # Amesos
      # AztecOO
      Belos
      # Epetra
      # Ifpack
      # MueLu
      Teuchos
      Tpetra
      Thyra
      # ShyLU
      ShyLU_DD
      Xpetra
      Zoltan
      Zoltan2)
  target_link_libraries(feat-kernel-solver PRIVATE Trilinos::all_selected_libs)

  #target_compile_options(feat-kernel-solver PRIVATE ${Trilinos_CXX_COMPILER_FLAGS})

  ## if trilinos has mumps available, we have to link against the fortran headers, which requires fortran to be an activated language of the project
  if(TARGET MUMPS::all_libs)
    set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "" FORCE)
    enable_language(Fortran)
    find_package(MPI COMPONENTS Fortran C CXX)
    # target_link_libraries(feat-kernel-solver PRIVATE MUMPS::all_libs MPI::MPI_CXX MPI::MPI_Fortran mpi_mpifh gfortran)
    target_link_libraries(feat-kernel-solver PRIVATE MUMPS::all_libs MPI::MPI_CXX MPI::MPI_Fortran)
    target_link_libraries(feat-kernel-solver PRIVATE ${Trilinos_Fortran_IMPLICIT_LINK_LIBRARIES})
  endif()
  target_link_libraries(feat-kernel-solver PRIVATE $<TARGET_NAME_IF_EXISTS:UMFPACK::all_libs>)
  target_link_libraries(feat-kernel-solver PRIVATE $<TARGET_NAME_IF_EXISTS:Scotch::all_libs>)
  target_link_libraries(feat-kernel-solver PRIVATE $<TARGET_NAME_IF_EXISTS:SCALAPACK::all_libs>)
  target_link_libraries(feat-kernel-solver PRIVATE $<TARGET_NAME_IF_EXISTS:LAPACK::all_libs>)
  target_link_libraries(feat-kernel-solver PRIVATE $<TARGET_NAME_IF_EXISTS:BLAS::all_libs>)
  # target_link_libraries(TRILINOS_TPL INTERFACE $<TARGET_NAME_IF_EXISTS:METIS::all_libs>)
  # target_link_libraries(feat-kernel-solver PUBLIC TRILINOS_TPL)
endif (FEAT_HAVE_TRILINOS)

if (FEAT_HAVE_UMFPACK)
  find_package(UMFPACK 6.3.5 REQUIRED)
  target_link_libraries(feat-kernel-solver PRIVATE SuiteSparse::UMFPACK)
endif (FEAT_HAVE_UMFPACK)


# list of solver tests
set(test_list
  amavanka-test
  basic_solver-test
  cudss-test
  direct_sparse_solver-test
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
