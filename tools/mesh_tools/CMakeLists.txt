# CMakeLists for directory "tools/mesh_tools"

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# tools to build
ADD_EXECUTABLE(mesh-indexer mesh_indexer.cpp)
TARGET_LINK_LIBRARIES(mesh-indexer feat)

ADD_EXECUTABLE(mesh-extruder mesh_extruder.cpp)
TARGET_LINK_LIBRARIES(mesh-extruder feat)

ADD_EXECUTABLE(mesh2eps mesh_to_eps.cpp)
TARGET_LINK_LIBRARIES(mesh2eps feat)

ADD_EXECUTABLE(mesh2vtk mesh_to_vtk.cpp)
TARGET_LINK_LIBRARIES(mesh2vtk feat)

ADD_EXECUTABLE(tri2mesh tri_to_mesh.cpp)
TARGET_LINK_LIBRARIES (tri2mesh feat)
