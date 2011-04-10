# Top-Level FEAST CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAST )

include (CTest)
ENABLE_TESTING()

set(MEMORYCHECK_COMMAND "valgrind")
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full")

# BUILD_ID mechanism
INCLUDE ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid.cmake )

if (FEAST_BACKENDS_CUDA)
  FIND_PACKAGE(CUDA)
endif(FEAST_BACKENDS_CUDA)


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
ADD_SUBDIRECTORY( applications )
ADD_SUBDIRECTORY( tools )
