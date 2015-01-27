# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

if (MPI_CXX_INCLUDE_PATH)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

#top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

ADD_EXECUTABLE(noscarc noscarc.cpp)

TARGET_LINK_LIBRARIES(noscarc feast ${MPI_CXX_LIBRARIES})
if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(noscarc PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()
if(MPI_CXX_LINK_FLAGS)
  set_target_properties(noscarc PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()
