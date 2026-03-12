# CMakeLists for directory "applications/gendie"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add all executables
if(FEAT_HAVE_CGAL)
  add_executable(steady-voxel-gendie steady_voxel_gendie.cpp)
  add_executable(steady-voxel-gendie-index32 steady_voxel_gendie.cpp)
  add_executable(material-viscosity-table material_viscosity_table.cpp)

  target_compile_definitions(steady-voxel-gendie-index32 PRIVATE FEAT_GENDIE_APP_INDEX32=1)

  # define Q1T app list
  set(app_list
    steady-voxel-gendie
    steady-voxel-gendie-index32
    material-viscosity-table
  )

  foreach (app ${app_list})
    target_link_libraries(${app} PRIVATE feat)
  endforeach (app)
endif(FEAT_HAVE_CGAL)
