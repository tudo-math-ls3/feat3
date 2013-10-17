# Top-Level FEAST CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAST )

include (CTest)
ENABLE_TESTING()
configure_file (${FEAST_SOURCE_DIR}/CTestCustom.cmake.in ${FEAST_BINARY_DIR}/CTestCustom.cmake COPONLY)
configure_file (${FEAST_SOURCE_DIR}/CTestCustom.cmake.in ${FEAST_BINARY_DIR}/CTestCustom.ctest COPONLY)

find_program(MEMORYCHECK_COMMAND NAMES valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full" CACHE STRING "valgrind options")

if (FEAST_BACKENDS_CUDA)
  set (CUDA_VERBOSE_BUILD ON)
  FIND_PACKAGE(CUDA)
  set (CUDA_VERBOSE_BUILD ON)
  set(CUDA_NVCC_FLAGS "-arch sm_13" CACHE STRING "Semi-colon delimit multiple arguments." FORCE)
  #TODO  use CUDA_HOST_COMPILER to get rid of hard wired gcc host compilation
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)
endif(FEAST_BACKENDS_CUDA)

if (FEAST_MPI)
  find_package(MPI REQUIRED)
endif (FEAST_MPI)


# store all variables describing the build in a custom header file
# to use later on in the real code
configure_file (${FEAST_SOURCE_DIR}/feast_config.hpp.in ${FEAST_BINARY_DIR}/feast_config.hpp)

# our custom target to build (not execute) all tests in all directories
ADD_CUSTOM_TARGET (tests)

# our custom target to wipe out all cmake generated files
INCLUDE ( ${FEAST_SOURCE_DIR}/cmake_modules/target_distclean.cmake )

# our custom target to build the documentation via doxygen
# TODO Redirect doxygen output to BINARY_DIR/doc
ADD_CUSTOM_TARGET ( doc doxygen ${FEAST_SOURCE_DIR}/doxy_in/doxygen_config WORKING_DIRECTORY ${FEAST_SOURCE_DIR} )
if ( NOT (${FEAST_SOURCE_DIR} STREQUAL ${FEAST_BINARY_DIR}) AND UNIX)
    ADD_CUSTOM_command( TARGET doc POST_BUILD
      COMMAND mkdir ${FEAST_BINARY_DIR}/doc
      COMMAND mv ${FEAST_SOURCE_DIR}/doc/html ${FEAST_BINARY_DIR}/doc)
endif ( NOT (${FEAST_SOURCE_DIR} STREQUAL ${FEAST_BINARY_DIR}) AND UNIX)


# directories to scan recursively
ADD_SUBDIRECTORY( test_system )
ADD_SUBDIRECTORY( kernel )
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( applications )
