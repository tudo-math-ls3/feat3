# CmakeLists for directory "applications/meshopt_refinement-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(meshopt_refinement-app meshopt_refinement-app)
TARGET_LINK_LIBRARIES(meshopt_refinement-app feat)
