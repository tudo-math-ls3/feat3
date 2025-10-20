# CMakeLists for directory "stupro"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add navsto-k-e-bfs-2d application (2D version)
add_executable(navsto-k-e-bfs-2d domain_control.cpp stokes_level.cpp stokes_solver.cpp check_point.cpp time_stepping.cpp vtk_writer.cpp navsto_k_e_bfs_2d.cpp)
target_link_libraries(navsto-k-e-bfs-2d feat)
if (MPI_C_COMPILE_FLAGS)
  set_target_properties(navsto-k-e-bfs-2d PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif ()
if (MPI_C_LINK_FLAGS)
  set_target_properties(navsto-k-e-bfs-2d PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif ()
