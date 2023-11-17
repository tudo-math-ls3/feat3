# CMakeLists for directory "tools"

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# add subdirectories with custom CMakeLists
ADD_SUBDIRECTORY( mesh_tools )

# add tools to build
ADD_EXECUTABLE(cub-list cub_list/cub_list.cpp)
TARGET_LINK_LIBRARIES(cub-list feat)

ADD_EXECUTABLE(matrix-cond matrix_cond/matrix_cond.cpp)
TARGET_LINK_LIBRARIES(matrix-cond feat)

ADD_EXECUTABLE(matrix-info matrix_info/matrix_info.cpp)
TARGET_LINK_LIBRARIES(matrix-info feat)

ADD_EXECUTABLE(csr2mtx io/csr_to_mtx.cpp)
TARGET_LINK_LIBRARIES (csr2mtx feat)

ADD_EXECUTABLE(mtx2csr io/mtx_to_csr.cpp)
TARGET_LINK_LIBRARIES (mtx2csr feat)

ADD_EXECUTABLE(dm2mtx io/dm_to_mtx.cpp)
TARGET_LINK_LIBRARIES (dm2mtx feat)

ADD_EXECUTABLE(mtx2dm io/mtx_to_dm.cpp)
TARGET_LINK_LIBRARIES (mtx2dm feat)

ADD_EXECUTABLE(dv2exp io/dv_to_exp.cpp)
TARGET_LINK_LIBRARIES (dv2exp feat)

ADD_EXECUTABLE(exp2dv io/exp_to_dv.cpp)
TARGET_LINK_LIBRARIES (exp2dv feat)

ADD_EXECUTABLE(dv2mtx io/dv_to_mtx.cpp)
TARGET_LINK_LIBRARIES (dv2mtx feat)

ADD_EXECUTABLE(mtx2dv io/mtx_to_dv.cpp)
TARGET_LINK_LIBRARIES (mtx2dv feat)

ADD_EXECUTABLE(csr2tga io/csr_to_tga.cpp)
TARGET_LINK_LIBRARIES (csr2tga feat)

ADD_EXECUTABLE(mtx2tga io/mtx_to_tga.cpp)
TARGET_LINK_LIBRARIES (mtx2tga feat)

ADD_EXECUTABLE(resort-matrix io/resort_matrix.cpp)
TARGET_LINK_LIBRARIES (resort-matrix feat)

ADD_EXECUTABLE(mpi-info mpi_info/mpi_info.cpp)
TARGET_LINK_LIBRARIES(mpi-info feat)

ADD_EXECUTABLE(checkpoint-info checkpoint_info/checkpoint_info.cpp)
TARGET_LINK_LIBRARIES(checkpoint-info feat)

ADD_EXECUTABLE(stokes2vtk stokes2vtk/stokes2vtk.cpp)
TARGET_LINK_LIBRARIES(stokes2vtk feat)

if (FEAT_HAVE_CUDA)
  CUDA_ADD_EXECUTABLE(cuda-device-info cuda_device_info/cuda_device_info.cpp)
  SET_PROPERTY(TARGET cuda-device-info APPEND_STRING PROPERTY COMPILE_FLAGS " -isystem${CUDA_INCLUDE_DIRS}")
  #TARGET_LINK_LIBRARIES (cuda-device-info)
endif(FEAT_HAVE_CUDA)
