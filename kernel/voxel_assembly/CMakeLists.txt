# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( arch )

# list of lafem source files
set(kernel-voxel-assembly-list
  empty_voxel.cpp
  )

add_library(feat-kernel-voxel-assembly STATIC ${kernel-voxel-assembly-list})
if (FEAT_COTIRE)
  cotire(feat-kernel-voxel-assembly)
endif (FEAT_COTIRE)

target_link_libraries(feat-kernel-voxel-assembly
PUBLIC
  feat-dependency
PRIVATE
  feat-kernel-voxel-assembly-arch feat-kernel-util)

if (FEAT_HAVE_CUDA)
  target_link_libraries(feat-kernel-voxel-assembly PRIVATE feat-kernel-util-cuda)
endif (FEAT_HAVE_CUDA)

# list of lafem tests
set(test_list
  voxel_assembly-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET voxel_assembly_tests CUDA)
