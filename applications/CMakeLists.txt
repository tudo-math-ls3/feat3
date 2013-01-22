# CmakeLists for directory "applications"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)


#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

# applications to build
ADD_EXECUTABLE(foundation_system_test-app foundation_system_test-app.cpp)
TARGET_LINK_LIBRARIES(foundation_system_test-app feast)

ADD_EXECUTABLE(foundation_comm_test-app foundation_comm_test-app.cpp)
TARGET_LINK_LIBRARIES(foundation_comm_test-app feast)

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
