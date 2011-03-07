# Top-Level FEAST CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAST )

ENABLE_TESTING()

INCLUDE ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid.cmake )


#enable DEBUG define if the option is set
OPTION (DEBUG_BUILD "Build with debug support.")
if (DEBUG_BUILD)
  add_definitions (-DDEBUG)
  #message ("Building with debug support...")
endif (DEBUG_BUILD)

# our custom target to build (not execute) all tests in all directories
ADD_CUSTOM_TARGET (tests)
#our custom target to execute all compiled tests
#ADD_CUSTOM_TARGET (check COMMAND ${FEAST_SOURCE_DIR}/test_system/run_tests ${FEAST_BINARY_DIR} )
#ADD_DEPENDENCIES( check tests )
# our custom target to wipe out all cmake generated files
INCLUDE ( ${FEAST_SOURCE_DIR}/cmake_modules/target_distclean.cmake )
ADD_DEPENDENCIES( distclean clean )
# our custom target to build the documentation via doxygen
# TODO Redirect doxygen output to BINARY_DIR/doc
ADD_CUSTOM_TARGET ( doc doxygen ${FEAST_SOURCE_DIR}/doc/doxygen_config WORKING_DIRECTORY ${FEAST_SOURCE_DIR} )
if ( NOT (${FEAST_SOURCE_DIR} STREQUAL ${FEAST_BINARY_DIR}) AND UNIX)
    ADD_CUSTOM_command( TARGET doc POST_BUILD
      COMMAND mkdir ${FEAST_BINARY_DIR}/doc
      COMMAND mv ${FEAST_SOURCE_DIR}/doc/html ${FEAST_BINARY_DIR}/doc)
endif ( NOT (${FEAST_SOURCE_DIR} STREQUAL ${FEAST_BINARY_DIR}) AND UNIX)

ADD_SUBDIRECTORY( kernel )
ADD_SUBDIRECTORY( test_system )
#ADD_SUBDIRECTORY( applications )
