# CMakeLists for directory "area51"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 3.9)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build

ADD_EXECUTABLE(basis-dump basis_dump.cpp)
TARGET_LINK_LIBRARIES(basis-dump feat)

ADD_EXECUTABLE(cgal-analyse-mesh cgal_analyse_mesh.cpp)
TARGET_LINK_LIBRARIES(cgal-analyse-mesh feat)

ADD_EXECUTABLE(cgal-analyse-mesh-parallel cgal_analyse_mesh_parallel.cpp)
TARGET_LINK_LIBRARIES(cgal-analyse-mesh-parallel feat)

ADD_EXECUTABLE(dbg-alg-dof-parti dbg_alg_dof_parti.cpp)
TARGET_LINK_LIBRARIES(dbg-alg-dof-parti feat)

ADD_EXECUTABLE(dbg-amavanka dbg_amavanka.cpp)
TARGET_LINK_LIBRARIES(dbg-amavanka feat)

ADD_EXECUTABLE(dbg-cubature dbg_cubature.cpp)
TARGET_LINK_LIBRARIES(dbg-cubature feat)

ADD_EXECUTABLE(dbg-dist-file-io dbg_dist_file_io.cpp)
TARGET_LINK_LIBRARIES(dbg-dist-file-io feat)

ADD_EXECUTABLE(dbg-isoparam-1 dbg_isoparam1.cpp)
TARGET_LINK_LIBRARIES(dbg-isoparam-1 feat)

ADD_EXECUTABLE(dbg-isoparam-2 dbg_isoparam2.cpp)
TARGET_LINK_LIBRARIES(dbg-isoparam-2 feat)

ADD_EXECUTABLE(dbg-isoparam-3 dbg_isoparam3.cpp)
TARGET_LINK_LIBRARIES(dbg-isoparam-3 feat)

ADD_EXECUTABLE(dbg-mesh-width dbg_mesh_width.cpp)
TARGET_LINK_LIBRARIES(dbg-mesh-width feat)

ADD_EXECUTABLE(dbg-meshperm dbg_meshperm.cpp)
TARGET_LINK_LIBRARIES(dbg-meshperm feat)

ADD_EXECUTABLE(dbg-nlopt dbg_nlopt.cpp)
TARGET_LINK_LIBRARIES(dbg-nlopt feat)

ADD_EXECUTABLE(dbg-pmdcdsc dbg_pmdcdsc.cpp)
TARGET_LINK_LIBRARIES(dbg-pmdcdsc feat)

ADD_EXECUTABLE(dbg-recursive-parti dbg_recursive_parti.cpp)
TARGET_LINK_LIBRARIES(dbg-recursive-parti feat)

ADD_EXECUTABLE(dbg-sphere-fem dbg_sphere_fem.cpp)
TARGET_LINK_LIBRARIES(dbg-sphere-fem feat)

ADD_EXECUTABLE(dbg-stokes-3field dbg_stokes_3field.cpp)
TARGET_LINK_LIBRARIES(dbg-stokes-3field feat)

ADD_EXECUTABLE(dbg-trace-1 dbg_trace_1.cpp)
TARGET_LINK_LIBRARIES(dbg-trace-1 feat)

ADD_EXECUTABLE(dbg-trace-2 dbg_trace_2.cpp)
TARGET_LINK_LIBRARIES(dbg-trace-2 feat)

ADD_EXECUTABLE(dbg-stabilize dbg_stabilize.cpp)
TARGET_LINK_LIBRARIES(dbg-stabilize feat)

ADD_EXECUTABLE(dbg-tuple-matrix dbg_tuple_matrix.cpp)
TARGET_LINK_LIBRARIES(dbg-tuple-matrix feat)

ADD_EXECUTABLE(dbg-vanka dbg_vanka.cpp)
TARGET_LINK_LIBRARIES(dbg-vanka feat)

if(FEAT_HAVE_UMFPACK)
   ADD_EXECUTABLE(feti-seq-usquare feti_seq_usquare.cpp)
   TARGET_LINK_LIBRARIES(feti-seq-usquare feat)
   if(FEAT_HAVE_MPI)
      ADD_EXECUTABLE(feti-para-usquare feti_para_usquare.cpp)
      TARGET_LINK_LIBRARIES(feti-para-usquare feat)

      ADD_EXECUTABLE(feti-multigrid-usquare feti_multigrid_usquare.cpp)
      TARGET_LINK_LIBRARIES(feti-multigrid-usquare feat)
   endif(FEAT_HAVE_MPI)
endif(FEAT_HAVE_UMFPACK)

######################### feti_multigrid_usquare
if(FEAT_HAVE_UMFPACK)
if (FEAT_HAVE_MPI)
  ADD_TEST(sleep304 sleep 2)
  SET_PROPERTY(TEST sleep304 PROPERTY LABELS "mpi,sleep")

  ADD_TEST(feti-multigrid-mpi-16 ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target feti-multigrid-usquare
    --build-nocmake
    --build-noclean
    --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16
    ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/feti-multigrid-usquare --test ${MPIEXEC_POSTFLAGS})
  SET_PROPERTY(TEST feti-multigrid-mpi-16 PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST feti-multigrid-mpi-16 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif(FEAT_HAVE_MPI)
endif(FEAT_HAVE_UMFPACK)
