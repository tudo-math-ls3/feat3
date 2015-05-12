# Top-Level FEAST CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAST )

site_name( FEAST_HOSTNAME )

include (CTest)
ENABLE_TESTING()

if (FEAST_BACKENDS_CUDA)
  FIND_PACKAGE(CUDA)
  set(CUDA_NVCC_FLAGS "-g -lineinfo -arch=${FEAST_CUDA_ARCH} -O3 --restrict --use_fast_math --compiler-options -Wall;-g")
  if (FEAST_CUDA_VERBOSE)
    set(CUDA_VERBOSE_BUILD ON)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -v -Xptxas=-v")
  endif (FEAST_CUDA_VERBOSE)
  if (FEAST_DEBUG_MODE)
    set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -G --compiler-options=-rdynamic")
  endif(FEAST_DEBUG_MODE)
  if (FEAST_CUDA_HOST_COMPILER)
    set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -ccbin " ${FEAST_CUDA_HOST_COMPILER})
  endif (FEAST_CUDA_HOST_COMPILER)
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)

  if (FEAST_CUDAMEMCHECK)
    find_program(CUDA_MEMCHECK_EXE cuda-memcheck)
    set(CUDA_MEMCHECK_EXE ${CUDA_MEMCHECK_EXE} --leak-check full)
  endif (FEAST_CUDAMEMCHECK)

  if (CUDA_VERSION_MAJOR GREATER 6)
    set(FEAST_HAVE_CUSOLVER ON)
  endif (CUDA_VERSION_MAJOR GREATER 6)
endif(FEAST_BACKENDS_CUDA)

if (FEAST_MPI)
  find_package(MPI REQUIRED)
  execute_process(COMMAND "${MPIEXEC}" --version
    WORKING_DIRECTORY ${FEAST_SOURCE_DIR}
    OUTPUT_VARIABLE CMAKE_MPI_VERSION
    ERROR_VARIABLE CMAKE_MPI_VERSION
    )
  string(REPLACE "\n" " " CMAKE_MPI_VERSION ${CMAKE_MPI_VERSION})
endif (FEAST_MPI)

if (FEAST_VALGRIND)
  find_program(VALGRIND_EXE valgrind)
  set (VALGRIND_EXE ${VALGRIND_EXE} --fullpath-after=${FEAST_SOURCE_DIR} --leak-check=full
    --track-origins=yes --partial-loads-ok=yes)
  #--keep-stacktraces=alloc-and-free
endif (FEAST_VALGRIND)


# our custom target to build (not execute) all tests in all directories
ADD_CUSTOM_TARGET (tests)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${FEAST_SOURCE_DIR}/doxy_in/doxygen_config.in ${FEAST_BINARY_DIR}/doxy_in/doxygen_config @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${FEAST_BINARY_DIR}/doxy_in/doxygen_config
    WORKING_DIRECTORY ${FEAST_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

if (CMAKE_CXX_COMPILER_ARG1)
  string(STRIP ${CMAKE_CXX_COMPILER_ARG1} CMAKE_CXX_COMPILER_ARG1_STRIPED)
  find_path(CMAKE_CXX_COMPILER_ARG1_PATH ${CMAKE_CXX_COMPILER_ARG1_STRIPED})
endif (CMAKE_CXX_COMPILER_ARG1)

# directories to scan recursively
ADD_SUBDIRECTORY( thirdparty )
ADD_SUBDIRECTORY( kernel )
ADD_SUBDIRECTORY( test_system )
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( applications )
ADD_SUBDIRECTORY( benchmarks )
ADD_SUBDIRECTORY( tutorials )

# store all variables describing the build in a custom header file
# to use later on in the real code
configure_file (${FEAST_SOURCE_DIR}/feast_config.hpp.in ${FEAST_BINARY_DIR}/feast_config.hpp)
