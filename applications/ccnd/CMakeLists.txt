# CMakeLists for directory "applications/ccnd"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add all executables
add_executable(ccnd-steady-dfg95-2d ccnd_steady_dfg95.cpp)
add_executable(ccnd-steady-dfg95-2d-q1t ccnd_steady_dfg95.cpp)
add_executable(ccnd-steady-dfg95-3d ccnd_steady_dfg95.cpp)
add_executable(ccnd-steady-dricav-2d ccnd_steady_dricav.cpp)
add_executable(ccnd-steady-taygre-2d ccnd_steady_taygre.cpp)
add_executable(ccnd-steady-ypipe-2d ccnd_steady_ypipe.cpp)
add_executable(ccnd-unsteady-dfg95-2d ccnd_unsteady_dfg95.cpp)
add_executable(ccnd-unsteady-dfg95-3d ccnd_unsteady_dfg95.cpp)

# define 2D app list
set(app_list_2d
  ccnd-steady-dfg95-2d
  ccnd-steady-dfg95-2d-q1t
  ccnd-steady-dricav-2d
  ccnd-steady-taygre-2d
  ccnd-steady-ypipe-2d
  ccnd-unsteady-dfg95-2d
)

# define 3D app list
set(app_list_3d
  ccnd-steady-dfg95-3d
  ccnd-unsteady-dfg95-3d
)

# define isoparametric Q2/P1dc app list
set(app_list_iso
  ccnd-steady-dfg95-2d
  ccnd-unsteady-dfg95-2d
  ccnd-steady-taygre-2d
  ccnd-steady-ypipe-2d
  ccnd-steady-dfg95-3d
  ccnd-unsteady-dfg95-3d
)

#define Q1T/P0 app list
set(app_list_q1t
  ccnd-steady-dfg95-2d-q1t
)

# list of all applications
set(app_list ${app_list_2d} ${app_list_3d})

# set dimension
foreach (app ${app_list_2d})
  target_compile_definitions(${app} PUBLIC FEAT_CCND_APP_DIM=2)
endforeach (app)
foreach (app ${app_list_3d})
  target_compile_definitions(${app} PUBLIC FEAT_CCND_APP_DIM=3)
endforeach (app)

# set isoparametric trafo for Q2/P1dc apps
foreach (app ${app_list_iso})
  target_compile_definitions(${app} PUBLIC FEAT_CCND_APP_ISOPARAM)
endforeach (app)

# set Q1T/P0 spaces
foreach (app ${app_list_q1t})
  target_compile_definitions(${app} PUBLIC FEAT_CCND_APP_Q1T_P0)
endforeach (app)

# set the other compile options for all applications
foreach (app ${app_list})
  target_link_libraries(${app} feat)
  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(${app} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(${app} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()
endforeach (app)

#######################################################################################################################

# quadmath versions of selected applications

if (FEAT_HAVE_QUADMATH)
  add_executable(ccnd-steady-dfg95-2d-qp ccnd_steady_dfg95.cpp)
  add_executable(ccnd-steady-dfg95-3d-qp ccnd_steady_dfg95.cpp)

  target_link_libraries(ccnd-steady-dfg95-2d-qp feat)
  target_link_libraries(ccnd-steady-dfg95-3d-qp feat)

  target_compile_definitions(ccnd-steady-dfg95-2d-qp PUBLIC FEAT_CCND_APP_DIM=2 FEAT_CCND_APP_QUADMATH FEAT_CCND_APP_ISOPARAM)
  target_compile_definitions(ccnd-steady-dfg95-3d-qp PUBLIC FEAT_CCND_APP_DIM=3 FEAT_CCND_APP_QUADMATH FEAT_CCND_APP_ISOPARAM)

  if (MPI_C_COMPILE_FLAGS)
    set_target_properties(ccnd-steady-dfg95-2d-qp PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  endif ()
  if (MPI_C_LINK_FLAGS)
    set_target_properties(ccnd-steady-dfg95-2d-qp PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif ()
endif (FEAT_HAVE_QUADMATH)
