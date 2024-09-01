# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of arch source files
if (NOT FEAT_HAVE_CUDA)
set(kernel-voxel-assembly-arch-list
  burgers_assembler.cpp
  defo_assembler.cpp
  poisson_assembler.cpp
)
else (NOT FEAT_HAVE_CUDA)
set(kernel-voxel-assembly-arch-list
  empty_arch.cpp
)
endif (NOT FEAT_HAVE_CUDA)

if (NOT FEAT_HAVE_OMP)
  message(STATUS "OMP not found, missing CPU threading capability of voxel assembly.")
endif ()

add_library(feat-kernel-voxel-assembly-arch STATIC ${kernel-voxel-assembly-arch-list})
target_link_libraries(feat-kernel-voxel-assembly-arch PUBLIC feat-dependency feat-kernel-util)

if (FEAT_COTIRE)
  cotire(feat-kernel-voxel-assembly-arch)
endif (FEAT_COTIRE)
target_link_libraries(feat-kernel-voxel-assembly-arch PRIVATE feat-kernel-util)

if (FEAT_HAVE_CUDA)
  add_library(feat-kernel-voxel-assembly-arch-cuda STATIC
    burgers_assembler.cu
    defo_assembler.cu
    poisson_assembler.cu
  )
  target_link_libraries(feat-kernel-voxel-assembly-arch-cuda
  PUBLIC
    feat-dependency
  PRIVATE
    feat-kernel-util-cuda)
  ## get openmp flags and convert (potential list with ; to list with ,)
  set(CUDA_TMP_FLAGS "")
  if(OpenMP_CXX_FOUND)
    string(REPLACE ";" "," CUDA_TMP_FLAGS "${OpenMP_CXX_FLAGS}")
  endif(OpenMP_CXX_FOUND)
  target_compile_options(feat-kernel-voxel-assembly-arch-cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:$<$<BOOL:FEAT_HAVE_OMP>:-Xcompiler=${CUDA_TMP_FLAGS}>>) #-dc for relocatable device code...
  target_compile_options(feat-kernel-voxel-assembly-arch-cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-suppress-stack-size-warning>)
  target_link_libraries(feat-kernel-voxel-assembly-arch
  PRIVATE
    feat-kernel-voxel-assembly-arch-cuda)
endif (FEAT_HAVE_CUDA)
