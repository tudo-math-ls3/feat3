# CMakeLists for directory "tools/mesh_tools"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# tools to build
add_executable(mesh-indexer mesh_indexer.cpp)
target_link_libraries(mesh-indexer feat)

add_executable(mesh-extruder mesh_extruder.cpp)
target_link_libraries(mesh-extruder feat)

add_executable(mesh-partitioner mesh_partitioner.cpp)
target_link_libraries(mesh-partitioner feat)

add_executable(mesh-validator mesh_validator.cpp)
target_link_libraries(mesh-validator feat)

add_executable(mesh2eps mesh_to_eps.cpp)
target_link_libraries(mesh2eps feat)

add_executable(mesh2svg mesh_to_svg.cpp)
target_link_libraries(mesh2svg feat)

add_executable(mesh2vtk mesh_to_vtk.cpp)
target_link_libraries(mesh2vtk feat)

add_executable(mesh2tri mesh_to_tri.cpp)
target_link_libraries (mesh2tri feat)

add_executable(tri2mesh tri_to_mesh.cpp)
target_link_libraries (tri2mesh feat)
