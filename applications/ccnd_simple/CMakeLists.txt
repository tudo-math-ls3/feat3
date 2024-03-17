# CMakeLists for directory "applications/ccnd_simple"

# this needs to be set again for each directory
cmake_minimum_required(VERSION 3.9)

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# we need the fparser library for all these apps
if (FEAT_HAVE_FPARSER)

  # add simple-steady-01 application (2D version, iso-parametric)
  add_executable(ccnd-simple-steady-01-2d domain.cpp stokes_level.cpp steady_solver.cpp vtk_writer.cpp app_steady_01.cpp)
  target_link_libraries(ccnd-simple-steady-01-2d feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-steady-01-2d PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_SIMPLE_ISOPARAM")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-steady-01-2d PROPERTIES COMPILE_FLAGS "-DFEAT_CCND_SIMPLE_ISOPARAM")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(ccnd-simple-steady-01-2d PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS} -DFEAT_CCND_SIMPLE_ISOPARAM")
  endif ()

  # add simple-steady-01 application (3D version, iso-parametric)
  add_executable(ccnd-simple-steady-01-3d domain.cpp stokes_level.cpp steady_solver.cpp vtk_writer.cpp app_steady_01.cpp)
  target_link_libraries(ccnd-simple-steady-01-3d feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-steady-01-3d PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_SIMPLE_3D -DFEAT_CCND_SIMPLE_ISOPARAM")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-steady-01-3d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_3D -DFEAT_CCND_SIMPLE_ISOPARAM")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(ccnd-simple-steady-01-3d PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()

  # add simple-unsteady-01 application (2D version, iso-parametric)
  add_executable(ccnd-simple-unsteady-01-2d domain.cpp stokes_level.cpp steady_solver.cpp unsteady_solver.cpp vtk_writer.cpp app_unsteady_01.cpp)
  target_link_libraries(ccnd-simple-unsteady-01-2d feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-2d PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_SIMPLE_ISOPARAM")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-2d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_ISOPARAM")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-2d PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()

  # add simple-unsteady-01 application (3D version, iso-parametric)
  add_executable(ccnd-simple-unsteady-01-3d domain.cpp stokes_level.cpp steady_solver.cpp unsteady_solver.cpp vtk_writer.cpp app_unsteady_01.cpp)
  target_link_libraries(ccnd-simple-unsteady-01-3d feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-3d PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -DFEAT_CCND_SIMPLE_3D -DFEAT_CCND_SIMPLE_ISOPARAM")
  else(MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-3d PROPERTIES COMPILE_FLAGS " -DFEAT_CCND_SIMPLE_3D -DFEAT_CCND_SIMPLE_ISOPARAM")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(ccnd-simple-unsteady-01-3d PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()

  #####################################################################################################################

  # TODO add tests

  #####################################################################################################################

  # add a custom target for all these apps
  add_custom_target(ccnd_simple_apps DEPENDS ccnd-simple-steady-01-2d ccnd-simple-steady-01-3d ccnd-simple-unsteady-01-2d ccnd-simple-unsteady-01-3d)

endif (FEAT_HAVE_FPARSER)
