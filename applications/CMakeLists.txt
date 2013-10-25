# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)


#top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# applications to build
ADD_EXECUTABLE(foundation_system_test-app foundation_system_test-app.cpp)
ADD_EXECUTABLE(foundation_comm_test-app foundation_comm_test-app.cpp)
ADD_EXECUTABLE(foundation_synch_test-app foundation_synch_test-app.cpp)
ADD_EXECUTABLE(scarc_test-app scarc_test-app.cpp)
ADD_EXECUTABLE(poisson_test-app poisson_test-app.cpp)
ADD_EXECUTABLE(proto-stokes proto_stokes.cpp)
TARGET_LINK_LIBRARIES(foundation_system_test-app feast)
TARGET_LINK_LIBRARIES(foundation_comm_test-app feast)
TARGET_LINK_LIBRARIES(foundation_synch_test-app feast)
TARGET_LINK_LIBRARIES(scarc_test-app feast)
TARGET_LINK_LIBRARIES(poisson_test-app feast)
TARGET_LINK_LIBRARIES(proto-stokes feast)

if (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_comm_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_comm_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "mpirun" -np 2 ${FEAST_BINARY_DIR}/applications/foundation_comm_test-app)
  SET_PROPERTY(TEST foundation_comm_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_comm_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_comm_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_comm_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_comm_test-app")
  SET_PROPERTY(TEST foundation_comm_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_comm_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (NOT FEAST_SERIAL_MODE)

if (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_system_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_system_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "mpirun" -np 4 ${FEAST_BINARY_DIR}/applications/foundation_system_test-app)
  SET_PROPERTY(TEST foundation_system_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_system_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_system_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_system_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-systemand "${FEAST_BINARY_DIR}/applications/foundation_system_test-app")
  SET_PROPERTY(TEST foundation_system_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_system_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (NOT FEAST_SERIAL_MODE)

if (NOT FEAST_SERIAL_MODE)
  ADD_TEST(scarc_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "mpirun" -np 2 ${FEAST_BINARY_DIR}/applications/scarc_test-app)
  SET_PROPERTY(TEST scarc_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST scarc_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (NOT FEAST_SERIAL_MODE)
  ADD_TEST(scarc_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target scarc_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/applications/scarc_test-app")
  SET_PROPERTY(TEST scarc_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST scarc_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (NOT FEAST_SERIAL_MODE)

if (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_synch_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_synch_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "mpirun" -np 2 ${FEAST_BINARY_DIR}/applications/foundation_synch_test-app)
  SET_PROPERTY(TEST foundation_synch_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST foundation_synch_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (NOT FEAST_SERIAL_MODE)
  ADD_TEST(foundation_synch_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target foundation_synch_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/applications/foundation_synch_test-app")
  SET_PROPERTY(TEST foundation_synch_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST foundation_synch_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (NOT FEAST_SERIAL_MODE)

if (NOT FEAST_SERIAL_MODE)
  ADD_TEST(poisson_test-app_mpi ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "mpirun" -np 4 ${FEAST_BINARY_DIR}/applications/poisson_test-app)
  SET_PROPERTY(TEST poisson_test-app_mpi PROPERTY LABELS "mpi")
  SET_PROPERTY(TEST poisson_test-app_mpi PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
else (NOT FEAST_SERIAL_MODE)
  ADD_TEST(poisson_test-app_serial ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target poisson_test-app
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/applications/poisson_test-app")
  SET_PROPERTY(TEST poisson_test-app_serial PROPERTY LABELS "serial")
  SET_PROPERTY(TEST poisson_test-app_serial PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
endif (NOT FEAST_SERIAL_MODE)

# add all tests to applications_tests
ADD_CUSTOM_TARGET(applications_tests DEPENDS scarc_test-app foundation_system_test-app)

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests applications_tests)
