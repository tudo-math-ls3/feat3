# Top-Level FEAT CMakeLists.txt

cmake_minimum_required(VERSION 3.18.0)
if (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 0)))
  cmake_policy(SET CMP0054 NEW)
endif (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 0)))

if (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 15)))
  cmake_policy(SET CMP0091 NEW) ##sets global MSVC runtime settings...
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
endif (((CMAKE_MAJOR_VERSION GREATER 3)) OR ((CMAKE_MAJOR_VERSION EQUAL 3) AND (CMAKE_MINOR_VERSION GREATER 15)))

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

project( FEAT )

#wrapper for sequential ctest execution
set(FEAT_CTEST_RUNNER "")

#if required, overwrite internal cxx flags
if(FEAT_DELETE_BUILD_FLAGS)
  message(STATUS "Delete internal build-config flags")
  set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING "Opt Cxx flag" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING "Debug Cxx flag" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "" CACHE STRING "RelDeb Cxx flag" FORCE)
  set(CMAKE_CXX_FLAGS_MINSIZEREL "" CACHE STRING "MinSize Cxx flag" FORCE)
endif(FEAT_DELETE_BUILD_FLAGS)

site_name( FEAT_HOSTNAME )

include (CTest)
enable_testing()

find_package(Git)
if (GIT_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
    WORKING_DIRECTORY ${FEAT_SOURCE_DIR}
    OUTPUT_VARIABLE FEAT_GIT_SHA1
    ERROR_VARIABLE FEAT_GIT_SHA1)
  if (FEAT_GIT_SHA1 MATCHES "fatal")
    set(FEAT_GIT_SHA1 "")
  endif (FEAT_GIT_SHA1 MATCHES "fatal")
  string(STRIP "${FEAT_GIT_SHA1}" FEAT_GIT_SHA1)
endif (GIT_FOUND)

if (FEAT_COTIRE)
  include(${FEAT_SOURCE_DIR}/build_system/cotire.cmake)
endif (FEAT_COTIRE)

if (FEAT_HAVE_MKL)
  set(BLA_VENDOR Intel10_64ilp_seq)
  find_package(BLAS REQUIRED)
  add_link_options(${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
endif (FEAT_HAVE_MKL)

if (FEAT_HAVE_CUDA)
  #TODO: We could also use Debug builds and so on... for now we simply erase these...
  set(CMAKE_CUDA_FLAGS_RELEASE "" CACHE STRING "Opt CUDA flag" FORCE)
  set(CMAKE_CUDA_FLAGS_DEBUG "" CACHE STRING "Debug CUDA flag" FORCE)
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "" CACHE STRING "RelDeb CUDA flag" FORCE)
  set(CMAKE_CUDA_FLAGS_MINSIZEREL "" CACHE STRING "MinSize CUDA flag" FORCE)

  set(INTERNAL_CUDA_FLAGS "-g -std=c++17 -O3 --restrict --compiler-options -Wall,-g,-Wno-deprecated-declarations")
  if (FEAT_CUDA_VERBOSE)
    set(CUDA_VERBOSE_BUILD ON)
    set(INTERNAL_CUDA_FLAGS "${INTERNAL_CUDA_FLAGS} -v -Xptxas=-v")
  endif (FEAT_CUDA_VERBOSE)
  if (FEAT_DEBUG_MODE)
    set (INTERNAL_CUDA_FLAGS "${INTERNAL_CUDA_FLAGS} -G --compiler-options=-rdynamic")
  else (FEAT_DEBUG_MODE)
    set (INTERNAL_CUDA_FLAGS "${INTERNAL_CUDA_FLAGS} -lineinfo")
  endif(FEAT_DEBUG_MODE)
  if (FEAT_CUDA_HOST_COMPILER)
    set (INTERNAL_CUDA_FLAGS "${INTERNAL_CUDA_FLAGS} -ccbin ${FEAT_CUDA_HOST_COMPILER}")
  endif (FEAT_CUDA_HOST_COMPILER)

  #Important: Set cuda flag before enable language call
  set(CMAKE_CUDA_FLAGS "${INTERNAL_CUDA_FLAGS}" CACHE STRING "Cuda compile flag")
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set (CUDA_PROPAGATE_HOST_FLAGS OFF)

  if (FEAT_CUDAMEMCHECK)
    find_program(CUDA_MEMCHECK_EXE compute-sanitizer)
    set(CUDA_MEMCHECK_EXE ${CUDA_MEMCHECK_EXE} --leak-check full)
  endif (FEAT_CUDAMEMCHECK)

  if (FEAT_HAVE_HALFMATH)
    include_directories("${CUDAToolkit_INCLUDE_DIRS}")
  endif (FEAT_HAVE_HALFMATH)

endif(FEAT_HAVE_CUDA)

if (FEAT_HAVE_MPI)
  if (FEAT_MPI_CXX)
    set(MPI_CXX_COMPILER ${FEAT_MPI_CXX})
    set(MPI_C_COMPILER ${FEAT_MPI_C})
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

  if (MPI_C_COMPILER_INCLUDE_DIRS)
    include_directories(BEFORE "${MPI_C_COMPILER_INCLUDE_DIRS}")
  endif()
  if (MPI_CXX_COMPILER_INCLUDE_DIRS)
    include_directories(BEFORE "${MPI_CXX_COMPILER_INCLUDE_DIRS}")
  endif ()

  add_definitions( -DOMPI_SKIP_MPICXX )

  #set ctest runner to single process mpiexec execution
  set(FEAT_CTEST_RUNNER ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS})
endif (FEAT_HAVE_MPI)

if ("${FEAT_COMPILER_ID}" STREQUAL "clang")
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  set(THREADS_PREFER_PTHREAD_FLAG TRUE)
  find_package(Threads REQUIRED)
  add_link_options(${CMAKE_THREAD_LIBS_INIT})
endif ()
find_package(OpenMP)
if (OpenMP_FOUND)
  option(FEAT_HAVE_OMP "FEAT OpenMP support available" ON)
  add_compile_options(${OpenMP_CXX_FLAGS})
  add_link_options(${OpenMP_CXX_FLAGS})
endif (OpenMP_FOUND)

if (FEAT_VALGRIND)
  find_program(VALGRIND_EXE valgrind)
  set (VALGRIND_EXE ${VALGRIND_EXE} -v --suppressions=${FEAT_SOURCE_DIR}/valgrind.supp --fullpath-after=${FEAT_SOURCE_DIR} --leak-check=full --partial-loads-ok=yes)
  #--track-origins=yes
  #--keep-stacktraces=alloc-and-free
endif (FEAT_VALGRIND)


# include directories of third-party libraries

if (FEAT_TRILINOS_PATH)
  find_package(Trilinos
      REQUIRED
      PATHS "${FEAT_TRILINOS_PATH}/"
      COMPONENTS
      Amesos
      AztecOO
      Belos
      Epetra
      Ifpack
      MueLu
      Teuchos
      Tpetra
      Thyra
      ShyLU
      Xpetra
      Zoltan
      Zoltan2)
    include_directories(AFTER ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
endif (FEAT_TRILINOS_PATH)


if (FEAT_HAVE_ALGLIB)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ALGLIB}/cpp/src")
endif(FEAT_HAVE_ALGLIB)

if (FEAT_HAVE_BOOST)
  if (FEAT_PREFER_EXTERNAL_LIBS)
    find_package(Boost REQUIRED)
    include_directories(BEFORE "${Boost_INCLUDE_DIRS}")
  else (FEAT_PREFER_EXTERNAL_LIBS)
    include_directories(BEFORE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_BOOST}")
  endif (FEAT_PREFER_EXTERNAL_LIBS)
endif (FEAT_HAVE_BOOST)

if (FEAT_HAVE_CGAL)
  if (NOT FEAT_HAVE_BOOST)
    find_package(Boost REQUIRED COMPONENTS thread)
    include_directories(BEFORE "${Boost_INCLUDE_DIRS}")
    if (NOT Boost_FOUND)
      message("boost not found on your system! you may add the keyword boost to the build-id-string in your configure statement to install it.")
    endif (NOT Boost_FOUND)
  endif (NOT FEAT_HAVE_BOOST)

  if (NOT BOOST_HAS_THREADS AND NOT OpenMP_FOUND)
    message(FATAL_ERROR "cgal only supported with multithreading capabilities. Please provide a valid boost or openMP application.")
  endif ()
  if (FEAT_PREFER_EXTERNAL_LIBS)
    find_package(CGAL 5.5 REQUIRED)
    include_directories("${CGAL_INCLUDE_DIRS}")
  else (FEAT_PREFER_EXTERNAL_LIBS)
    include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_CGAL}/include")
  endif (FEAT_PREFER_EXTERNAL_LIBS)
endif(FEAT_HAVE_CGAL)

if (FEAT_HAVE_DEATH_HANDLER)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_DEATHHANDLER}/DeathHandler-master")
endif(FEAT_HAVE_DEATH_HANDLER)

if (FEAT_HAVE_FLOATX)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_FLOATX}/FloatX-develop/src")
endif(FEAT_HAVE_FLOATX)

if (FEAT_HAVE_FPARSER)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_FPARSER}/")
endif(FEAT_HAVE_FPARSER)

if (FEAT_HAVE_PARMETIS)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/include")
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/include")
endif (FEAT_HAVE_PARMETIS)

if (FEAT_HAVE_TRIANGLE)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_TRIANGLE}")
endif(FEAT_HAVE_TRIANGLE)

if (FEAT_HAVE_ZFP)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZFP}/include/")
endif(FEAT_HAVE_ZFP)

if (FEAT_HAVE_ZLIB)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZLIB}/")
endif(FEAT_HAVE_ZLIB)

if (FEAT_HAVE_ZOLTAN)
  include_directories("${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/include")
endif(FEAT_HAVE_ZOLTAN)

# our custom target to build (but not execute) all tests in all directories
add_custom_target(tests)

# add a target to generate API documentation with Doxygen
if (${CMAKE_VERSION} VERSION_GREATER "3.8.99")
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    set(DOXYGEN_ABBREVIATE_BRIEF "The $name class" "The $name widget" "The $name file" is provides specifies contains represents a an the)
    set(DOXYGEN_ALIASES
      [["compilerhack=\xrefitem compilerhacks \"Compiler Hack\" \"Compiler Hacks\""]]
      [["platformswitch=\xrefitem platformswitches \"Platform dependend code branching\" \"Platform depenendend code branches\""]]
      [["resident=\ref res_vs_tran_resident \"resident\""]]
      [["transient=\ref res_vs_tran_transient \"transient\""]])

    set(DOXYGEN_CITE_BIB_FILES "${FEAT_SOURCE_DIR}/doxy_in/references.bib")
    set(DOXYGEN_COLLABORATION_GRAPH NO)
    set(DOXYGEN_EXCLUDE_PATTERNS "*/applications/*" "*/area51/*" "*/benchmarks/*" "*/build_system/*" "*/data/*" "*/thirdparty/*" "*/tools/*" "*/tutorials/*" "*-test*.cpp")
    set(DOXYGEN_EXPAND_AS_DEFINED "DOXY" "static_assert")
    set(DOXYGEN_EXTRA_PACKAGES "amsmath" "amssymb")
    set(DOXYGEN_FILE_PATTERS "*.c" "*.cc" "*.cxx" "*.cpp" "*.c++" "*.cu" "*.dox" "*.py" "*.h" "*.hh" "*.hpp" "*.h++" "*.f90" "*.f")
    set(DOXYGEN_FORMULA_FONTSIZE 12)
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
include_directories("${FEAT_BINARY_DIR}" "${FEAT_SOURCE_DIR}")

# directories to scan recursively
add_subdirectory( thirdparty )
add_subdirectory( kernel )
add_subdirectory( test_system )
add_subdirectory( tools )
add_subdirectory( applications )
add_subdirectory( area51 )
add_subdirectory( benchmarks )
add_subdirectory( tutorials )
add_subdirectory( control )

# store all variables describing the build in a custom header file
# to use later on in the real code
configure_file(${FEAT_SOURCE_DIR}/feat_config.hpp.in ${FEAT_BINARY_DIR}/feat_config.hpp)
