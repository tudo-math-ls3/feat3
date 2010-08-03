# Top-Level FEAST CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAST )

# our custom target to build (not execute) all tests in all directories
ADD_CUSTOM_TARGET (tests)
# our custom target to wipe out all cmake generated files
INCLUDE( ${FEAST_SOURCE_DIR}/cmake_modules/target_distclean.cmake OPTIONAL )

ADD_SUBDIRECTORY( kernel )
#ADD_SUBDIRECTORY( test_system )
#ADD_SUBDIRECTORY( applications )
