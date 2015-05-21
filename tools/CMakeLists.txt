# CmakeLists for directory "tools"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})


# applications to build
ADD_EXECUTABLE(basis-dump basis_dump/basis_dump.cpp)
TARGET_LINK_LIBRARIES(basis-dump feast)

ADD_EXECUTABLE(cub-list cub_list/cub_list.cpp)
TARGET_LINK_LIBRARIES(cub-list feast)

ADD_EXECUTABLE(dbg-bcasm dbg_bcasm/dbg_bcasm.cpp)
TARGET_LINK_LIBRARIES(dbg-bcasm feast)

ADD_EXECUTABLE(dbg-fe-asm-1 dbg_fe_asm_1/dbg_fe_asm_1.cpp)
TARGET_LINK_LIBRARIES(dbg-fe-asm-1 feast)

ADD_EXECUTABLE(dbg-fe-asm-2 dbg_fe_asm_2/dbg_fe_asm_2.cpp)
TARGET_LINK_LIBRARIES(dbg-fe-asm-2 feast)

ADD_EXECUTABLE(dbg-mirror dbg_mirror/dbg_mirror.cpp)
TARGET_LINK_LIBRARIES(dbg-mirror feast)

ADD_EXECUTABLE(dbg-nodefunc dbg_nodefunc/dbg_nodefunc.cpp)
TARGET_LINK_LIBRARIES(dbg-nodefunc feast)

ADD_EXECUTABLE(dbg-patch dbg_patch/dbg_patch.cpp)
TARGET_LINK_LIBRARIES(dbg-patch feast)

ADD_EXECUTABLE(matrix-info matrix_info/matrix_info.cpp)
TARGET_LINK_LIBRARIES(matrix-info feast)

ADD_EXECUTABLE(csr2mtx io/csr_to_mtx.cpp)
TARGET_LINK_LIBRARIES (csr2mtx feast)

ADD_EXECUTABLE(ell2mtx io/ell_to_mtx.cpp)
TARGET_LINK_LIBRARIES (ell2mtx feast)

ADD_EXECUTABLE(mesh_streamer-app mesh_streamer/mesh_streamer-app.cpp)
TARGET_LINK_LIBRARIES (mesh_streamer-app feast)

ADD_EXECUTABLE(mtx2csr io/mtx_to_csr.cpp)
TARGET_LINK_LIBRARIES (mtx2csr feast)

ADD_EXECUTABLE(mtx2ell io/mtx_to_ell.cpp)
TARGET_LINK_LIBRARIES (mtx2ell feast)

ADD_EXECUTABLE(csr2ell io/csr_to_ell.cpp)
TARGET_LINK_LIBRARIES (csr2ell feast)

ADD_EXECUTABLE(ell2csr io/ell_to_csr.cpp)
TARGET_LINK_LIBRARIES (ell2csr feast)

ADD_EXECUTABLE(dv2exp io/dv_to_exp.cpp)
TARGET_LINK_LIBRARIES (dv2exp feast)

ADD_EXECUTABLE(exp2dv io/exp_to_dv.cpp)
TARGET_LINK_LIBRARIES (exp2dv feast)

ADD_EXECUTABLE(csr2tga io/csr_to_tga.cpp)
TARGET_LINK_LIBRARIES (csr2tga feast)

ADD_EXECUTABLE(ell2tga io/ell_to_tga.cpp)
TARGET_LINK_LIBRARIES (ell2tga feast)

ADD_EXECUTABLE(mtx2tga io/mtx_to_tga.cpp)
TARGET_LINK_LIBRARIES (mtx2tga feast)

ADD_EXECUTABLE(resort-csr io/resort_csr.cpp)
TARGET_LINK_LIBRARIES (resort-csr feast)
