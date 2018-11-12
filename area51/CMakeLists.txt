# CMakeLists for directory "area51"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build

ADD_EXECUTABLE(basis-dump basis_dump.cpp)
TARGET_LINK_LIBRARIES(basis-dump feat)

ADD_EXECUTABLE(dbg-alg-dof-parti dbg_alg_dof_parti.cpp)
TARGET_LINK_LIBRARIES(dbg-alg-dof-parti feat)

ADD_EXECUTABLE(dbg-amavanka dbg_amavanka.cpp)
TARGET_LINK_LIBRARIES(dbg-amavanka feat)

ADD_EXECUTABLE(dbg-mesh-width dbg_mesh_width.cpp)
TARGET_LINK_LIBRARIES(dbg-mesh-width feat)

ADD_EXECUTABLE(dbg-nlopt dbg_nlopt.cpp)
TARGET_LINK_LIBRARIES(dbg-nlopt feat)

ADD_EXECUTABLE(dbg-pmdcdsc dbg_pmdcdsc.cpp)
TARGET_LINK_LIBRARIES(dbg-pmdcdsc feat)

ADD_EXECUTABLE(dbg-recursive-parti dbg_recursive_parti.cpp)
TARGET_LINK_LIBRARIES(dbg-recursive-parti feat)

ADD_EXECUTABLE(dbg-stokes-3field dbg_stokes_3field.cpp)
TARGET_LINK_LIBRARIES(dbg-stokes-3field feat)

ADD_EXECUTABLE(dbg-trace-1 dbg_trace_1.cpp)
TARGET_LINK_LIBRARIES(dbg-trace-1 feat)

ADD_EXECUTABLE(dbg-trace-2 dbg_trace_2.cpp)
TARGET_LINK_LIBRARIES(dbg-trace-2 feat)

ADD_EXECUTABLE(dbg-streamdiff dbg_streamdiff.cpp)
TARGET_LINK_LIBRARIES(dbg-streamdiff feat)

ADD_EXECUTABLE(dbg-tuple-matrix dbg_tuple_matrix.cpp)
TARGET_LINK_LIBRARIES(dbg-tuple-matrix feat)

ADD_EXECUTABLE(dbg-vanka dbg_vanka.cpp)
TARGET_LINK_LIBRARIES(dbg-vanka feat)
