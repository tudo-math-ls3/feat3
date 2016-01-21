# CmakeLists for directory "applications/meshopt_boundary-app"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(meshopt_boundary-app meshopt_boundary-app)
TARGET_LINK_LIBRARIES(meshopt_boundary-app feast)
