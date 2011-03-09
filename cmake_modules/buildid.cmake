#add build_id logic to the build system

set (BUILD_ID "" CACHE STRING "Use specific build id or NONE or empty for autodetection")
if (BUILD_ID STREQUAL "NONE")
  #do nothing, use env variable or whatever the user provides
elseif (BUILD_ID STREQUAL "")
  #use guess id
  #until then we enforce basic options
  #set (CMAKE_CXX_COMPILER "mpic++")
  #set (CMAKE_CXX_COMPILER_ARG1 "")
  #set (CMAKE_CXX_FLAGS "-O0")
  #set (CMAKE_SHARED_LINKER_FLAGS "")
  #set (CMAKE_EXE_LINKER_FLAGS "")
  #set (CMAKE_MODULE_LINKER_FLAGS "")
  #set (DEBUG_BUILD ON)
else ()
  #use specific id
  message ("Build ID not known!")
endif ()
