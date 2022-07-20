# CmakeLists for directory "tools"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 3.9)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(cub-list cub_list/cub_list.cpp)
TARGET_LINK_LIBRARIES(cub-list feat)

ADD_EXECUTABLE(matrix-cond matrix_cond/matrix_cond.cpp)
TARGET_LINK_LIBRARIES(matrix-cond feat)

ADD_EXECUTABLE(matrix-info matrix_info/matrix_info.cpp)
TARGET_LINK_LIBRARIES(matrix-info feat)

ADD_EXECUTABLE(mesh-indexer mesh_indexer/mesh_indexer.cpp)
TARGET_LINK_LIBRARIES(mesh-indexer feat)

ADD_EXECUTABLE(mesh-extruder mesh_extruder/mesh_extruder.cpp)
TARGET_LINK_LIBRARIES(mesh-extruder feat)

ADD_EXECUTABLE(mesh2eps mesh_to_eps/mesh_to_eps.cpp)
TARGET_LINK_LIBRARIES(mesh2eps feat)

ADD_EXECUTABLE(mesh2vtk mesh_to_vtk/mesh_to_vtk.cpp)
TARGET_LINK_LIBRARIES(mesh2vtk feat)

ADD_EXECUTABLE(csr2mtx io/csr_to_mtx.cpp)
TARGET_LINK_LIBRARIES (csr2mtx feat)

ADD_EXECUTABLE(ell2mtx io/ell_to_mtx.cpp)
TARGET_LINK_LIBRARIES (ell2mtx feat)

ADD_EXECUTABLE(mtx2csr io/mtx_to_csr.cpp)
TARGET_LINK_LIBRARIES (mtx2csr feat)

ADD_EXECUTABLE(mtx2ell io/mtx_to_ell.cpp)
TARGET_LINK_LIBRARIES (mtx2ell feat)

ADD_EXECUTABLE(csr2ell io/csr_to_ell.cpp)
TARGET_LINK_LIBRARIES (csr2ell feat)

ADD_EXECUTABLE(ell2csr io/ell_to_csr.cpp)
TARGET_LINK_LIBRARIES (ell2csr feat)

ADD_EXECUTABLE(dv2exp io/dv_to_exp.cpp)
TARGET_LINK_LIBRARIES (dv2exp feat)

ADD_EXECUTABLE(exp2dv io/exp_to_dv.cpp)
TARGET_LINK_LIBRARIES (exp2dv feat)

ADD_EXECUTABLE(csr2tga io/csr_to_tga.cpp)
TARGET_LINK_LIBRARIES (csr2tga feat)

ADD_EXECUTABLE(ell2tga io/ell_to_tga.cpp)
TARGET_LINK_LIBRARIES (ell2tga feat)

ADD_EXECUTABLE(mtx2tga io/mtx_to_tga.cpp)
TARGET_LINK_LIBRARIES (mtx2tga feat)

ADD_EXECUTABLE(resort-matrix io/resort_matrix.cpp)
TARGET_LINK_LIBRARIES (resort-matrix feat)

ADD_EXECUTABLE(tri2mesh io/tri_to_mesh.cpp)
TARGET_LINK_LIBRARIES (tri2mesh feat)

ADD_EXECUTABLE(mpi-info mpi_info/mpi_info.cpp)
TARGET_LINK_LIBRARIES(mpi-info feat)

ADD_EXECUTABLE(checkpoint-info checkpoint_info/checkpoint_info.cpp)
TARGET_LINK_LIBRARIES(checkpoint-info feat)

if (FEAT_HAVE_CUDA)
  CUDA_ADD_EXECUTABLE(cuda-device-info cuda_device_info/cuda_device_info.cpp)
  SET_PROPERTY(TARGET cuda-device-info APPEND_STRING PROPERTY COMPILE_FLAGS " -isystem${CUDA_INCLUDE_DIRS}")
  TARGET_LINK_LIBRARIES (cuda-device-info)
endif(FEAT_HAVE_CUDA)

if (FEAT_HAVE_FPARSER)
  ADD_EXECUTABLE(mesh-parsed-hit-test-vtk mesh_parsed_hit_test_vtk/mesh_parsed_hit_test_vtk.cpp)
  TARGET_LINK_LIBRARIES(mesh-parsed-hit-test-vtk feat)
endif(FEAT_HAVE_FPARSER)
