# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# include test_aux subdirectory
add_subdirectory( test_aux )

# list of geometry source files and depenencies
set(kernel-geometry-list
  parti_parmetis.cpp
  parti_zoltan.cpp
  voxel_map.cpp
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(cgal.cpp PROPERTIES COMPILE_OPTIONS "-Wno-maybe-uninitialized")
endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

if (FEAT_EICKT)
  set(kernel-geometry-list
    ${kernel-geometry-list}
    conformal_mesh-eickt.cpp
    mesh_part-eickt.cpp
    mesh_node-eickt.cpp
    mesh_file_reader-eickt.cpp
  )
endif (FEAT_EICKT)


add_library(feat-kernel-geometry STATIC ${kernel-geometry-list})

if(FEAT_HAVE_PARMETIS)
  find_package(parmetis 4.0.3 REQUIRED)
  target_link_libraries(feat-kernel-geometry PRIVATE ParMETIS::ParMETIS)
endif(FEAT_HAVE_PARMETIS)

if(FEAT_HAVE_ZOLTAN)
  find_package(zoltan 3.901 REQUIRED)
  target_link_libraries(feat-kernel-geometry PRIVATE Zoltan::Zoltan)
endif()

if(FEAT_HAVE_FPARSER)
  find_package(fparser 4.5.2 REQUIRED)
  target_link_libraries(feat-kernel-geometry PRIVATE fparser::fparser)
endif()

# CGAL wrapper library to hide CGALs dependency on boost
# We build the cgal-wrapper as a static library linking to the external cgal
# dependency and parts of FEAT. This static library is then linked privately
# into the geometry kernel. That way cgals dependencies on boost are hidden.
if(FEAT_HAVE_CGAL)
  find_package(CGAL 5.5.5 REQUIRED)
  add_library(feat-kernel-geometry-cgal STATIC cgal.cpp)

  target_link_libraries(feat-kernel-geometry-cgal
  PRIVATE
    CGAL::CGAL
  PUBLIC # TODO: necessary?
    feat-dependency feat-kernel-util)

  target_link_libraries(feat-kernel-geometry PRIVATE feat-kernel-geometry-cgal)
endif()

if (FEAT_VALGRIND)
  target_compile_definitions(feat-kernel-geometry PRIVATE CGAL_DISABLE_ROUNDING_MATH_CHECK=ON)
endif (FEAT_VALGRIND)

target_link_libraries(feat-kernel-geometry
PUBLIC
  feat-dependency
PRIVATE
  feat-kernel-util)

# list of geometry tests
set(test_list
  adaptive_refinement-test
  boundary_factory-test
  cgal-test
  hit_test_factory-test
  index_calculator-test
  mesh_node-test-conf-quad
  mesh_part-test
  shape_convert-test
  standard_refinery-test-conf-quad
  standard_refinery-test-conf-hexa
  standard_refinery-test-conf-tria
  standard_refinery-test-conf-tetra
  standard_refinery-test-struct-quad
  standard_refinery-test-struct-hexa
  struct_index_mapping-test
  parsed_hit_test_factory-test
  template_builder-test
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET geometry_tests LINK_AGAINST "feat-kernel-geometry-test-aux")

#disabling warnings with pragmas does not always work with gcc
#https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(cgal.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-undef")
  set_source_files_properties(cgal-test.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-undef")

  #silence weird glibcxx warning
  #cgal/CGAL-4.14/include/CGAL/AABB_tree.h:1138:20:
  #warning: argument 1 value ‘18446744073709551615’ exceeds maximum object size 9223372036854775807 [-Walloc-size-larger-than=]
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0)
    set_source_files_properties(cgal.cpp APPEND_STRING PROPERTY COMPILE_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    set_source_files_properties(cgal-test.cpp APPEND_STRING PROPERTY COMPILE_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    set_source_files_properties(cgal.cpp APPEND_STRING PROPERTY LINK_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
    set_source_files_properties(cgal-test.cpp APPEND_STRING PROPERTY LINK_OPTIONS " -Wno-alloc-size-larger-than -Wno-undef")
  endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0)
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
