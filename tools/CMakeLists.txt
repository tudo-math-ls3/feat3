# CMakeLists for directory "tools"

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# add subdirectories with custom CMakeLists
add_subdirectory( mesh_tools )

# add tools to build
add_executable(cub-list cub_list/cub_list.cpp)
target_link_libraries(cub-list feat)

add_executable(matrix-cond matrix_cond/matrix_cond.cpp)
target_link_libraries(matrix-cond feat)

add_executable(matrix-info matrix_info/matrix_info.cpp)
target_link_libraries(matrix-info feat)

add_executable(csr2mtx io/csr_to_mtx.cpp)
target_link_libraries(csr2mtx feat)

add_executable(mtx2csr io/mtx_to_csr.cpp)
target_link_libraries(mtx2csr feat)

add_executable(dm2mtx io/dm_to_mtx.cpp)
target_link_libraries(dm2mtx feat)

add_executable(mtx2dm io/mtx_to_dm.cpp)
target_link_libraries(mtx2dm feat)

add_executable(dv2exp io/dv_to_exp.cpp)
target_link_libraries(dv2exp feat)

add_executable(exp2dv io/exp_to_dv.cpp)
target_link_libraries(exp2dv feat)

add_executable(dv2mtx io/dv_to_mtx.cpp)
target_link_libraries(dv2mtx feat)

add_executable(mtx2dv io/mtx_to_dv.cpp)
target_link_libraries(mtx2dv feat)

add_executable(csr2tga io/csr_to_tga.cpp)
target_link_libraries(csr2tga feat)

add_executable(mtx2tga io/mtx_to_tga.cpp)
target_link_libraries(mtx2tga feat)

add_executable(resort-matrix io/resort_matrix.cpp)
target_link_libraries(resort-matrix feat)

add_executable(mpi-info mpi_info/mpi_info.cpp)
target_link_libraries(mpi-info feat)

add_executable(checkpoint-info checkpoint_info/checkpoint_info.cpp)
target_link_libraries(checkpoint-info feat)

add_executable(stokes2vtk stokes2vtk/stokes2vtk.cpp)
target_link_libraries(stokes2vtk feat)

if (FEAT_HAVE_CUDA)
  add_executable(cuda-device-info cuda_device_info/cuda_device_info.cpp)
  set_property(TARGET cuda-device-info APPEND_STRING PROPERTY COMPILE_FLAGS " -isystem${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
  target_link_libraries(cuda-device-info ${CUDA_cudart_LIBRARY})
endif (FEAT_HAVE_CUDA)
