# Top-Level FEAT CMakeLists.txt

cmake_minimum_required (VERSION 3.9)
if (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION LESS 13)))
  cmake_policy(SET CMP0058 OLD)
endif (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION LESS 13)))

if (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 0)))
cmake_policy(SET CMP0054 NEW)
endif (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 0)))

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

PROJECT( FEAT )

site_name( FEAT_HOSTNAME )

include (CTest)
ENABLE_TESTING()

find_package(Git)
if(GIT_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
    WORKING_DIRECTORY ${FEAT_SOURCE_DIR}
    OUTPUT_VARIABLE FEAT_GIT_SHA1
    ERROR_VARIABLE FEAT_GIT_SHA1)
  if (FEAT_GIT_SHA1 MATCHES "fatal")
    SET(FEAT_GIT_SHA1 "")
  endif (FEAT_GIT_SHA1 MATCHES "fatal")
  STRING(STRIP "${FEAT_GIT_SHA1}" FEAT_GIT_SHA1)
endif(GIT_FOUND)

if (FEAT_COTIRE)
  include(${FEAT_SOURCE_DIR}/build_system/cotire.cmake)
endif (FEAT_COTIRE)

if (FEAT_HAVE_CUDA)
  FIND_PACKAGE(CUDA)
  set(CUDA_NVCC_FLAGS " -g -arch=${FEAT_CUDA_ARCH} -std=c++11 -O3 --restrict --compiler-options -Wall,-g")
  if(BUILD_ID MATCHES ".*fast.*")
    set(CUDA_NVCC_FLAGS " --use_fast_math ${CUDA_NVCC_FLAGS}")
  endif(BUILD_ID MATCHES ".*fast.*")
  if (FEAT_CUDA_VERBOSE)
    set(CUDA_VERBOSE_BUILD ON)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -v -Xptxas=-v")
  endif (FEAT_CUDA_VERBOSE)
  if (FEAT_DEBUG_MODE)
    set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -G --compiler-options=-rdynamic")
  else (FEAT_DEBUG_MODE)
    set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -lineinfo")
  endif(FEAT_DEBUG_MODE)
  if (FEAT_CUDA_HOST_COMPILER)
    set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -ccbin " ${FEAT_CUDA_HOST_COMPILER})
  endif (FEAT_CUDA_HOST_COMPILER)
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)

  if (FEAT_CUDAMEMCHECK)
    find_program(CUDA_MEMCHECK_EXE cuda-memcheck)
    set(CUDA_MEMCHECK_EXE ${CUDA_MEMCHECK_EXE} --leak-check full)
  endif (FEAT_CUDAMEMCHECK)

  if (CMAKE_MAJOR_VERSION GREATER 2)
    if (CUDA_cusolver_LIBRARY)
      set(FEAT_HAVE_CUSOLVER ON)
    endif (CUDA_cusolver_LIBRARY)
  else (CMAKE_MAJOR_VERSION GREATER 2)
    if (CUDA_VERSION_MAJOR GREATER 6)
      set(FEAT_HAVE_CUSOLVER ON)
    endif (CUDA_VERSION_MAJOR GREATER 6)
  endif (CMAKE_MAJOR_VERSION GREATER 2)

  if (CUDA_VERSION_MAJOR GREATER 7)
    #set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " --compiler-options=-fopenmp")
  endif (CUDA_VERSION_MAJOR GREATER 7)
endif(FEAT_HAVE_CUDA)

if (FEAT_HAVE_MPI)
  if (FEAT_MPI_CXX)
    set (MPI_CXX_COMPILER ${FEAT_MPI_CXX})
    set (MPI_C_COMPILER ${FEAT_MPI_C})
  endif (FEAT_MPI_CXX)
  find_package(MPI REQUIRED)
  execute_process(COMMAND "${MPIEXEC}" --version
    WORKING_DIRECTORY ${FEAT_SOURCE_DIR}
    OUTPUT_VARIABLE CMAKE_MPI_VERSION
    ERROR_VARIABLE CMAKE_MPI_VERSION
    )
  if (CMAKE_MPI_VERSION)
    string(REPLACE "\n" " " CMAKE_MPI_VERSION ${CMAKE_MPI_VERSION})
  endif (CMAKE_MPI_VERSION)

  if (MPI_C_INCLUDE_PATH)
    include_directories(${MPI_C_INCLUDE_PATH})
  endif()

  ADD_DEFINITIONS( -DOMPI_SKIP_MPICXX )
endif (FEAT_HAVE_MPI)

find_package(OpenMP)
if (OpenMP_FOUND)
  OPTION(FEAT_HAVE_OMP "FEAT OpenMP support available" ON)
  ADD_COMPILE_OPTIONS(${OpenMP_CXX_FLAGS})
  ADD_LINK_OPTIONS(${OpenMP_CXX_FLAGS})
endif (OpenMP_FOUND)

if (FEAT_VALGRIND)
  find_program(VALGRIND_EXE valgrind)
  set (VALGRIND_EXE ${VALGRIND_EXE} -v --suppressions=${FEAT_SOURCE_DIR}/valgrind.supp --fullpath-after=${FEAT_SOURCE_DIR} --leak-check=full --partial-loads-ok=yes)
  #--track-origins=yes
  #--keep-stacktraces=alloc-and-free
endif (FEAT_VALGRIND)


# include directories of third-party libraries

if (FEAT_HAVE_ALGLIB)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ALGLIB}/cpp/src")
endif(FEAT_HAVE_ALGLIB)

if (FEAT_HAVE_CGAL)
  #set(Boost_USE_STATIC_LIBS ON)
  find_package(Boost REQUIRED COMPONENTS thread)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_CGAL}/include")
endif(FEAT_HAVE_CGAL)

if (FEAT_HAVE_DEATH_HANDLER)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_DEATHHANDLER}/DeathHandler-master")
endif(FEAT_HAVE_DEATH_HANDLER)

if (FEAT_HAVE_FLOATX)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_FLOATX}/FloatX-develop/src")
endif(FEAT_HAVE_FLOATX)

if (FEAT_HAVE_FPARSER)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_FPARSER}/")
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_HALFMATH)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HALF}/include")
endif(FEAT_HAVE_HALFMATH)

if (FEAT_HAVE_PARMETIS)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/include")
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/include")
endif (FEAT_HAVE_PARMETIS)

if (FEAT_HAVE_TRIANGLE)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_TRIANGLE}")
endif(FEAT_HAVE_TRIANGLE)

if (FEAT_HAVE_ZFP)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZFP}/include/")
endif(FEAT_HAVE_ZFP)

if (FEAT_HAVE_ZLIB)
  INCLUDE_DIRECTORIES("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZLIB}/")
endif(FEAT_HAVE_ZLIB)

# our custom target to build (but not execute) all tests in all directories
ADD_CUSTOM_TARGET (tests)

# add a target to generate API documentation with Doxygen
if (${CMAKE_VERSION} VERSION_GREATER "3.8.99")
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    set(DOXYGEN_ABBREVIATE_BRIEF "The $name class" "The $name widget" "The $name file" is provides specifies contains represents a an the)
    set(DOXYGEN_ALIASES
      [["compilerhack=\xrefitem compilerhacks \"Compiler Hack\" \"Compiler Hacks\""]]
      [["platformswitch=\xrefitem platformswitches \"Platform dependend code branching\" \"Platform depenendend code branches\""]])
    set(DOXYGEN_CITE_BIB_FILES "${FEAT_SOURCE_DIR}/doxy_in/references.bib")
    set(DOXYGEN_COLLABORATION_GRAPH NO)
    set(DOXYGEN_EXCLUDE_PATTERNS "*/applications/*" "*/area51/*" "*/benchmarks/*" "*/build_system/*" "*/data/*" "*/thirdparty/*" "*/tools/*" "*/tutorials/*" "*-test*.cpp")
    set(DOXYGEN_EXPAND_AS_DEFINED "DOXY" "static_assert")
    set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb")
    set(DOXYGEN_FILE_PATTERS "*.c" "*.cc" "*.cxx" "*.cpp" "*.c++" "*.cu" "*.dox" "*.py" "*.h" "*.hh" "*.hpp" "*.h++" "*.f90" "*.f")
    set(DOXYGEN_FORMULA_FONTSIZE 11)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_MAN NO)
    set(DOXYGEN_HTML_TIMESTAMP YES)
    set(DOXYGEN_IMAGE_PATH "${FEAT_SOURCE_DIR}/doxy_in")
    set(DOXYGEN_INLINE_INHERITED_MEMB YES)
    set(DOXYGEN_OUTPUT_DIRECTORY "${FEAT_BINARY_DIR}/doc")
    set(DOXYGEN_PREDEFINED "DOXYGEN=1" "DOXY(x)=x" "static_assert(c,m)=")
    set(DOXYGEN_TAB_SIZE "2")
    set(DOXYGEN_VERBATIM_VARS DOXYGEN_ALIASES)
    set(DOXYGEN_WARN_LOGFILE "${FEAT_BINARY_DIR}/doc/warnings.log")
    doxygen_add_docs(doc ${FEAT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen")
    set(DOXYGEN_ENABLED_SECTIONS "internal")
    doxygen_add_docs(doc-internal ${FEAT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen")
  endif (DOXYGEN_FOUND)
endif (${CMAKE_VERSION} VERSION_GREATER "3.8.99")

if (CMAKE_CXX_COMPILER_ARG1)
  string(STRIP ${CMAKE_CXX_COMPILER_ARG1} CMAKE_CXX_COMPILER_ARG1_STRIPED)
  find_path(CMAKE_CXX_COMPILER_ARG1_PATH ${CMAKE_CXX_COMPILER_ARG1_STRIPED})
endif (CMAKE_CXX_COMPILER_ARG1)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAT_BINARY_DIR}" "${FEAT_SOURCE_DIR}")

# directories to scan recursively
ADD_SUBDIRECTORY( thirdparty )
ADD_SUBDIRECTORY( kernel )
ADD_SUBDIRECTORY( test_system )
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( applications )
ADD_SUBDIRECTORY( area51 )
ADD_SUBDIRECTORY( benchmarks )
ADD_SUBDIRECTORY( tutorials )
ADD_SUBDIRECTORY( control )

# store all variables describing the build in a custom header file
# to use later on in the real code
configure_file (${FEAT_SOURCE_DIR}/feat_config.hpp.in ${FEAT_BINARY_DIR}/feat_config.hpp)
