cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

if (FEAST_HAVE_ALGLIB)
  # list of lafem source files
  file(GLOB alglib-list "${FEAST_SOURCE_DIR}/thirdparty/alglib/src/*.cpp")

  ADD_LIBRARY(thirdparty-alglib STATIC ${alglib-list})
endif (FEAST_HAVE_ALGLIB)
