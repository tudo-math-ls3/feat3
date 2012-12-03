# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)


#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

# applications to build
ADD_EXECUTABLE(dbg-cub dbg_cub/dbg_cub.cpp)
TARGET_LINK_LIBRARIES(dbg-cub feast)

ADD_EXECUTABLE(dbg-fe-asm-1 dbg_fe_asm_1/dbg_fe_asm_1.cpp)
TARGET_LINK_LIBRARIES(dbg-fe-asm-1 feast)

ADD_EXECUTABLE(dbg-mirror dbg_mirror/dbg_mirror.cpp)
TARGET_LINK_LIBRARIES(dbg-mirror feast)

ADD_EXECUTABLE(dbg-bcasm dbg_bcasm/dbg_bcasm.cpp)
TARGET_LINK_LIBRARIES(dbg-bcasm feast)

ADD_EXECUTABLE(dbg-nodefunc dbg_nodefunc/dbg_nodefunc.cpp)
TARGET_LINK_LIBRARIES(dbg-nodefunc feast)
