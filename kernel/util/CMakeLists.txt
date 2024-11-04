# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

set(KERNEL_LIBRARY_FILES
  dist.cpp
  dist_file_io.cpp
  kahan_summation.cpp
  memory_pool.cpp
  property_map.cpp
  statistics.cpp
  xml_scanner.cpp
)

if (WIN32)
  set(KERNEL_LIBRARY_FILES ${KERNEL_LIBRARY_FILES} os_windows.cpp)
endif (WIN32)

# list of util source files
add_library(kernel-util STATIC ${KERNEL_LIBRARY_FILES})

if (WIN32)
  target_link_libraries(kernel-util PRIVATE -limagehlp)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_link_libraries(kernel-util PRIVATE -lmsvcr120)
  endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
endif (WIN32)

# if(DEFINED ENV{OS})
#   #check if we are on Windows
#   if("$ENV{OS}" STREQUAL "Windows_NT")
#     #we probably only need to do this with libc based linkers, i.e. check if we use gcc
#     if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
#       TARGET_LINK_LIBRARIES(kernel-util -limagehlp)
#       TARGET_LINK_LIBRARIES(kernel-util -lmsvcr120)
#     endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
#     if(${MINGW})
#     if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
#       TARGET_LINK_LIBRARIES(kernel-util -limagehlp)
#       TARGET_LINK_LIBRARIES(kernel-util -lmsvcr120)
#     endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
#     endif(${MINGW})
#   endif("$ENV{OS}" STREQUAL "Windows_NT")
# endif(DEFINED ENV{OS})

if (FEAT_COTIRE)
  cotire(kernel-util)
endif (FEAT_COTIRE)

#link against mpi if available
target_link_libraries(kernel-util PRIVATE ${MPI_C_LIBRARIES})

if (FEAT_USE_LIKWID)
  if (FEAT_EXPLICIT_LIKWID_PATH)
    message (STATUS "Found likwid in ${FEAT_EXPLICIT_LIKWID_PATH}")
    set (likwid_FOUND ON)
    set (LIKWID_LIBRARIES ${FEAT_EXPLICIT_LIKWID_PATH}/lib/liblikwid.so)
    set (LIKWID_INCLUDES ${FEAT_EXPLICIT_LIKWID_PATH}/include)
  else (FEAT_EXPLICIT_LIKWID_PATH)
    find_package(likwid)
  endif (FEAT_EXPLICIT_LIKWID_PATH)
  if (likwid_FOUND)
    # set(LIKWID_LIBRARIES -L/sfw/likwid/git230818/lib -llikwid)
    target_link_libraries(kernel-util PUBLIC ${LIKWID_LIBRARIES})
    target_include_directories(kernel-util PUBLIC ${LIKWID_INCLUDES})
    # target_link_options(kernel-util PUBLIC LIKWID_PERFMON) # set by our config script
  else (likwid_FOUND)
    message(STATUS "Could not find likwid cmake path, consider specifing it with likwid_DIR")
  endif (likwid_FOUND)
endif (FEAT_USE_LIKWID)

if (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
  target_link_libraries(kernel-util PRIVATE -lexecinfo)
endif ()

if (MPI_C_COMPILE_FLAGS)
  set_target_properties(kernel-util PROPERTIES
    COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif ()
if (MPI_C_LINK_FLAGS)
  set_target_properties(kernel-util PROPERTIES
    LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif ()

if (FEAT_HAVE_CUDA)
  add_library(kernel-util-cuda STATIC
    cuda_util.cu
  )
  target_link_libraries(kernel-util-cuda PRIVATE ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY})
  target_link_libraries(kernel-util PRIVATE kernel-util-cuda)
endif (FEAT_HAVE_CUDA)

# list of util tests
set(test_list
  binary_stream-test
  math-test
  memory_usage-test
  meta_math-test
  pack-test
  property_map-test
  random-test
  simple_arg_parser-test
  string-test
  string_mapped-test
  tiny_algebra-test
  xml_scanner-test
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} feat test_system)

  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/util/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/util/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)
endforeach (test)

# add all tests to util_tests
add_custom_target(util_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests util_tests)
