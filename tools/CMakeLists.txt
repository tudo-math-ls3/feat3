# CmakeLists for directory "tools"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(basis-dump basis_dump/basis_dump.cpp)
TARGET_LINK_LIBRARIES(basis-dump feast)

ADD_EXECUTABLE(cub-list cub_list/cub_list.cpp)
TARGET_LINK_LIBRARIES(cub-list feast)

ADD_EXECUTABLE(dbg-fe-asm-1 dbg_fe_asm_1/dbg_fe_asm_1.cpp)
TARGET_LINK_LIBRARIES(dbg-fe-asm-1 feast)

ADD_EXECUTABLE(dbg-fe-asm-2 dbg_fe_asm_2/dbg_fe_asm_2.cpp)
TARGET_LINK_LIBRARIES(dbg-fe-asm-2 feast)

ADD_EXECUTABLE(dbg-nlopt dbg_nlopt/dbg_nlopt.cpp)
TARGET_LINK_LIBRARIES(dbg-nlopt feast)

ADD_EXECUTABLE(dbg-vanka dbg_vanka/dbg_vanka.cpp)
TARGET_LINK_LIBRARIES(dbg-vanka feast)

ADD_EXECUTABLE(matrix-cond matrix_cond/matrix_cond.cpp)
TARGET_LINK_LIBRARIES(matrix-cond feast)

ADD_EXECUTABLE(matrix-info matrix_info/matrix_info.cpp)
TARGET_LINK_LIBRARIES(matrix-info feast)

ADD_EXECUTABLE(mesh-extruder mesh_extruder/mesh_extruder.cpp)
TARGET_LINK_LIBRARIES(mesh-extruder feast)

ADD_EXECUTABLE(mesh-preprocessor mesh_preprocessor/mesh_preprocessor.cpp)
TARGET_LINK_LIBRARIES(mesh-preprocessor feast)

ADD_EXECUTABLE(mesh2vtk mesh_to_vtk/mesh_to_vtk.cpp)
TARGET_LINK_LIBRARIES(mesh2vtk feast)

ADD_EXECUTABLE(csr2mtx io/csr_to_mtx.cpp)
TARGET_LINK_LIBRARIES (csr2mtx feast)

ADD_EXECUTABLE(ell2mtx io/ell_to_mtx.cpp)
TARGET_LINK_LIBRARIES (ell2mtx feast)

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

ADD_EXECUTABLE(resort-matrix io/resort_matrix.cpp)
TARGET_LINK_LIBRARIES (resort-matrix feast)
