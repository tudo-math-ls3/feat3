# CmakeLists for directory "benchmarks"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})


# applications to build
ADD_EXECUTABLE(product-matvec product_matvec.cpp)
TARGET_LINK_LIBRARIES(product-matvec feast)
