# CMakeLists for directory "applications/ccnd_simple"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add simple-steady-01 application (2D version, iso-parametric)
add_executable(ccnd-simple-steady-01-2d domain_control.cpp stokes_level.cpp stokes_solver.cpp vtk_writer.cpp app_steady_01.cpp)
target_link_libraries(ccnd-simple-steady-01-2d feat)
set_target_properties(ccnd-simple-steady-01-2d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_DIM=2 -DFEAT_CCND_SIMPLE_ISOPARAM")

# add simple-steady-01 application (3D version, iso-parametric)
add_executable(ccnd-simple-steady-01-3d domain_control.cpp stokes_level.cpp stokes_solver.cpp vtk_writer.cpp app_steady_01.cpp)
target_link_libraries(ccnd-simple-steady-01-3d feat)
set_target_properties(ccnd-simple-steady-01-3d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_DIM=3 -DFEAT_CCND_SIMPLE_ISOPARAM")

# add simple-unsteady-01 application (2D version, iso-parametric)
add_executable(ccnd-simple-unsteady-01-2d domain_control.cpp stokes_level.cpp stokes_solver.cpp check_point.cpp time_stepping.cpp vtk_writer.cpp app_unsteady_01.cpp)
target_link_libraries(ccnd-simple-unsteady-01-2d feat)
set_target_properties(ccnd-simple-unsteady-01-2d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_DIM=2 -DFEAT_CCND_SIMPLE_ISOPARAM")

# add simple-unsteady-01 application (3D version, iso-parametric)
add_executable(ccnd-simple-unsteady-01-3d domain_control.cpp stokes_level.cpp stokes_solver.cpp check_point.cpp time_stepping.cpp vtk_writer.cpp app_unsteady_01.cpp)
target_link_libraries(ccnd-simple-unsteady-01-3d feat)
set_target_properties(ccnd-simple-unsteady-01-3d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_DIM=3 -DFEAT_CCND_SIMPLE_ISOPARAM")

#####################################################################################################################

# add a custom target for all these apps
add_custom_target(ccnd-simple-apps DEPENDS ccnd-simple-steady-01-2d ccnd-simple-steady-01-3d ccnd-simple-unsteady-01-2d ccnd-simple-unsteady-01-3d)
