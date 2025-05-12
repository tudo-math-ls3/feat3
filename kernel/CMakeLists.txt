# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# include kernel subdirectories
add_subdirectory( adjacency )
add_subdirectory( analytic )
add_subdirectory( assembly )
add_subdirectory( cubature )
add_subdirectory( geometry )
add_subdirectory( global )
add_subdirectory( lafem )
add_subdirectory( meshopt )
add_subdirectory( space )
add_subdirectory( solver )
add_subdirectory( trafo )
add_subdirectory( util )
add_subdirectory( voxel_assembly )

# list of kernel source files
set(KERNEL_LIBRARY_FILES
  backend.cpp
  runtime.cpp
)

# add kernel-root library
add_library(feat-kernel-root STATIC ${KERNEL_LIBRARY_FILES})
target_link_libraries(feat-kernel-root PUBLIC feat-dependency)
if (FEAT_FINALIZE_RESETS_DEVICE)
  target_compile_definitions(feat-kernel-root PRIVATE FEAT_FINALIZE_RESETS_DEVICE)
endif (FEAT_FINALIZE_RESETS_DEVICE)

if (FEAT_TESTING_VC)
  target_compile_definitions(feat-kernel-root PRIVATE FEAT_TESTING_VC)
endif (FEAT_TESTING_VC)

# feat-kernel-root and kernel-util have cyclic dependency
target_link_libraries(feat-kernel-root PRIVATE feat-kernel-util)
target_link_libraries(feat-kernel-util PRIVATE feat-kernel-root)

# add main feat library
add_library(feat INTERFACE)

# add additional alias to be cmake conformant
add_library(feat::feat ALIAS feat)

# list of libraries to include in feat lib
target_link_libraries(feat INTERFACE
  feat-kernel-adjacency
# feat-kernel-assembly
  feat-kernel-cubature
  feat-kernel-geometry
  feat-kernel-lafem
  feat-kernel-meshopt
# feat-kernel-space
  feat-kernel-solver
# feat-kernel-trafo
  feat-kernel-util
  feat-kernel-voxel-assembly
  feat-kernel-root
)

if (FEAT_HAVE_DEATH_HANDLER)
  find_package(death_handler REQUIRED)
  target_link_libraries(feat-kernel-root PRIVATE deathhandler::deathhandler)
endif (FEAT_HAVE_DEATH_HANDLER)

if (FEAT_HAVE_FPARSER)
  find_package(fparser 4.5.2 REQUIRED)
  target_link_libraries(feat INTERFACE fparser::fparser)
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_TRILINOS)
  #LINK_DIRECTORIES( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
  ## TODO: clear this dependency
  target_link_libraries(feat INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif (FEAT_HAVE_TRILINOS)

#libquadmath needs to be the last library linked against
if (FEAT_HAVE_QUADMATH)
  target_link_libraries(feat INTERFACE -lquadmath)
endif (FEAT_HAVE_QUADMATH)

if(FEAT_HAVE_PMP)
  find_package(pmp 3.0.0 REQUIRED)
  target_link_libraries(feat-kernel-root PUBLIC pmp::pmp)
endif()

# list of kernel tests
set(test_list
  shape-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET kernel_tests)
