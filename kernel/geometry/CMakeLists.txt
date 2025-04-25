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
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  add_dependencies(${test} feat-kernel-geometry-test-aux)
  target_link_libraries(${test} PRIVATE feat test_system feat-kernel-geometry-test-aux)


  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/geometry/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/geometry/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

endforeach (test)


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

# add all tests to geometry_tests
add_custom_target(geometry_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests geometry_tests)
