# CmakeLists for directory "applications/navier_stokes_ppnd"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(navier_stokes_ppnd navier_stokes_ppnd)
TARGET_LINK_LIBRARIES(navier_stokes_ppnd feat)
