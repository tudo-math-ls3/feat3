# CMakeLists for directory "area51"

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( ccnd_fiber )


# applications to build

add_executable(dbg-dist-mesh-distortion dbg_dist_mesh_distortion.cpp)
target_link_libraries(dbg-dist-mesh-distortion feat)

add_executable(dbg-mesh-distortion dbg_mesh_distortion.cpp)
target_link_libraries(dbg-mesh-distortion feat)

add_executable(dbg-voxel-domain dbg_voxel_domain.cpp)
target_link_libraries(dbg-voxel-domain feat)

add_executable(dbg-voxel-domain-gendie dbg_voxel_domain_gendie.cpp)
target_link_libraries(dbg-voxel-domain-gendie feat)

add_executable(basis-dump basis_dump.cpp)
target_link_libraries(basis-dump feat)

add_executable(cgal-analyse-mesh cgal_analyse_mesh.cpp)
target_link_libraries(cgal-analyse-mesh feat)

add_executable(cgal-analyse-mesh-parallel cgal_analyse_mesh_parallel.cpp)
target_link_libraries(cgal-analyse-mesh-parallel feat)

add_executable(dbg-alg-dof-parti dbg_alg_dof_parti.cpp)
target_link_libraries(dbg-alg-dof-parti feat)

add_executable(dbg-amavanka dbg_amavanka.cpp)
target_link_libraries(dbg-amavanka feat)

add_executable(dbg-base-splitter dbg_base_splitter.cpp)
target_link_libraries(dbg-base-splitter feat)

add_executable(dbg-cubature dbg_cubature.cpp)
target_link_libraries(dbg-cubature feat)

add_executable(dbg-dist-file-io dbg_dist_file_io.cpp)
target_link_libraries(dbg-dist-file-io feat)

add_executable(dbg-isoparam-1 dbg_isoparam1.cpp)
target_link_libraries(dbg-isoparam-1 feat)

add_executable(dbg-isoparam-2 dbg_isoparam2.cpp)
target_link_libraries(dbg-isoparam-2 feat)

add_executable(dbg-isoparam-3 dbg_isoparam3.cpp)
target_link_libraries(dbg-isoparam-3 feat)

add_executable(dbg-integrator dbg_integrator.cpp)
target_link_libraries(dbg-integrator feat)

add_executable(dbg-mesh-width dbg_mesh_width.cpp)
target_link_libraries(dbg-mesh-width feat)

add_executable(dbg-meshperm dbg_meshperm.cpp)
target_link_libraries(dbg-meshperm feat)

add_executable(dbg-nlopt dbg_nlopt.cpp)
target_link_libraries(dbg-nlopt feat)

add_executable(dbg-pmdcdsc dbg_pmdcdsc.cpp)
target_link_libraries(dbg-pmdcdsc feat)

add_executable(dbg-recursive-parti dbg_recursive_parti.cpp)
target_link_libraries(dbg-recursive-parti feat)

add_executable(dbg-slip-asm dbg_slip_asm.cpp)
target_link_libraries(dbg-slip-asm feat)

add_executable(dbg-sphere-fem dbg_sphere_fem.cpp)
target_link_libraries(dbg-sphere-fem feat)

add_executable(dbg-stokes-3field dbg_stokes_3field.cpp)
target_link_libraries(dbg-stokes-3field feat)

add_executable(dbg-trace-1 dbg_trace_1.cpp)
target_link_libraries(dbg-trace-1 feat)

add_executable(dbg-trace-2 dbg_trace_2.cpp)
target_link_libraries(dbg-trace-2 feat)

add_executable(dbg-trace-3 dbg_trace_3.cpp)
target_link_libraries(dbg-trace-3 feat)

add_executable(dbg-stabilize dbg_stabilize.cpp)
target_link_libraries(dbg-stabilize feat)

add_executable(dbg-tuple-matrix dbg_tuple_matrix.cpp)
target_link_libraries(dbg-tuple-matrix feat)

add_executable(dbg-vanka dbg_vanka.cpp)
target_link_libraries(dbg-vanka feat)

add_executable(dbg-parti-domain-control dbg_parti_domain_control.cpp)
target_link_libraries(dbg-parti-domain-control feat)

add_executable(hfem-direct hfem_direct.cpp)
target_link_libraries(hfem-direct feat)

if (FEAT_HAVE_UMFPACK)
   add_executable(feti-seq-usquare feti_seq_usquare.cpp)
   target_link_libraries(feti-seq-usquare feat)
   if (FEAT_HAVE_MPI)
      add_executable(feti-para-usquare feti_para_usquare.cpp)
      target_link_libraries(feti-para-usquare feat)

      add_executable(feti-multigrid-usquare feti_multigrid_usquare.cpp)
      target_link_libraries(feti-multigrid-usquare feat)
   endif (FEAT_HAVE_MPI)
endif (FEAT_HAVE_UMFPACK)


if (FEAT_HAVE_UMFPACK)
  if (FEAT_HAVE_MPI)
    add_test(sleep304 sleep 2)
    set_property(TEST sleep304 PROPERTY LABELS "mpi,sleep")

    add_test(feti-multigrid-mpi-16 ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target feti-multigrid-usquare
      --build-nocmake
      --build-noclean
      --test-command ${MPIEXEC} --map-by node ${MPIEXEC_NUMPROC_FLAG} 16
      ${MPIEXEC_PREFLAGS} ${FEAT_BINARY_DIR}/area51/feti-multigrid-usquare --test ${MPIEXEC_POSTFLAGS})
    set_property(TEST feti-multigrid-mpi-16 PROPERTY LABELS "mpi")
    set_property(TEST feti-multigrid-mpi-16 PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_HAVE_MPI)
endif (FEAT_HAVE_UMFPACK)
