# define minimum required cmake
# Set to 3.28 for EXCLUDE_FROM_ALL support in FetchContent_declare
cmake_minimum_required(VERSION 3.28.0)

# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

if(NOT CMAKE_BUILD_TYPE)
  if(FEAT_DEBUG_MODE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "FEAT3 build type")
  else()
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "FEAT3 build type")
  endif()
endif()

# our project
project(
  FEAT
  VERSION 1.0
  DESCRIPTION "Finite Element Analysis Toolbox"
  LANGUAGES C CXX)

# Add our own modules to cmake's module path
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/build_system/cmake)

#wrapper for sequential ctest execution
set(FEAT_CTEST_RUNNER "")

# if we are the main project (i.e. not included by add_subdirectories), we now do a few setup things
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # flag identifiying us as main project
  set(FEAT_MAIN_PROJECT YES INTERNAL)

  # set cxx standard we use across the board
  set(CMAKE_CXX_STANDARD 17)

  # Let's nicely support folders in IDEs <- TODO have to test if this break anything?
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # also, as main project, we have to include testing
  include(CTest)
  # no call to enable_testing() reuqired...

  # do git things
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

  # cotire only sensible if we are main project?
  option(FEAT_COTIRE "Add cotire compilation support to FEAT" OFF)
  if (FEAT_COTIRE)
    include(${FEAT_SOURCE_DIR}/build_system/cotire.cmake)
  endif (FEAT_COTIRE)

  # also, if main project, enable documentation
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
    set(DOXYGEN_ENABLE_PREPROCESSING YES)
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
    set(DOXYGEN_MACRO_EXPANSION YES)
    set(DOXYGEN_OUTPUT_DIRECTORY "${FEAT_BINARY_DIR}/doc")
    set(DOXYGEN_PREDEFINED "DOXYGEN=1" "DOXY(x)=x" "static_assert(c,m)=")
    set(DOXYGEN_TAB_SIZE "2")
    set(DOXYGEN_VERBATIM_VARS DOXYGEN_ALIASES)
    set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
    set(DOXYGEN_WARN_LOGFILE "${FEAT_BINARY_DIR}/doc/warnings.log")
    doxygen_add_docs(doc ${FEAT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen")
    set(DOXYGEN_ENABLED_SECTIONS "internal")
    doxygen_add_docs(doc-internal ${FEAT_SOURCE_DIR} COMMENT "Generating API documentation with Doxygen")
  endif (DOXYGEN_FOUND)
endif()

#INTERFACE target to add general build dependencies
add_library(feat-dependency INTERFACE EXCLUDE_FROM_ALL)

target_include_directories(feat-dependency INTERFACE ${FEAT_SOURCE_DIR} ${FEAT_BINARY_DIR})

# are we in debug mode? this should be handled in another way, if we want to have multiple build directories, i.e. IDE support
option(FEAT_DEBUG_MODE "Debug modus?" OFF)

# find and enable core dependencies
option(FEAT_HAVE_MKL "Support MKL backend for FEAT?" OFF)
if(FEAT_HAVE_MKL)
  # TODO: More general blas support
  set(BLA_VENDOR Intel10_64ilp_seq)
  find_package(BLAS REQUIRED)
  target_link_libraries(feat-dependency INTERFACE BLAS::BLAS) #TODO: necessary here? Should it not be enough
  #to link the actual implementations???
endif(FEAT_HAVE_MKL)

# cuda support?
option(FEAT_HAVE_CUDA "Enable cuda language and cuda backend?" OFF)
if(FEAT_HAVE_CUDA)
  # since we call this before enable_language, this should be set as long as the user does not provide them
  set(CMAKE_CUDA_FLAGS "-std=c++17 --restrict --compiler-options -Wall,-g,-Wno-deprecated-declarations" CACHE STRING "General CUDA FLAGS")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -O3" CACHE STRING "Cuda rel with deb")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3" CACHE STRING "Cuda rel with deb")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -O1" CACHE STRING "Cuda rel with deb")

  option(FEAT_CUDA_VERBOSE "Enable verbose cuda compiler?" OFF)
  if(FEAT_CUDA_VERBOSE)
    set(CUDA_VERBOSE_BUILD ON)
    target_compile_options(feat-dependency INTERFACE $<$<COMPILE_LANGUAGE:CUDA>:-v -Xptxas=-v>)
  endif(FEAT_CUDA_VERBOSE)

  if(FEAT_DEBUG_MODUS)
    target_compile_options(feat-dependency INTERFACE $<$<COMPILE_LANGUAGE:CUDA>:-G --compiler-options=-rdynamic>)
  else(FEAT_DEBUG_MODUS)
    target_compile_options(feat-dependency INTERFACE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)
  endif(FEAT_DEBUG_MODUS)

  ## host compiler <- thios should be done another way, correct, since this can simply be done
  # by setting CMAKE_CUDA_HOST_COMPILER

  option(FEAT_CUDA_ALLOW_UNSUPPORTED_COMPILER "Allow non officially tested compilers?" OFF)
  if(FEAT_CUDA_ALLOW_UNSUPPORTED_COMPILER)
    target_compile_options(feat-dependency INTERFACE $<$<COMPILE_LANGUAGE:CUDA>:--allow-unsupported-compiler>)
  endif(FEAT_CUDA_ALLOW_UNSUPPORTED_COMPILER)


  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)

  option(FEAT_CUDAMEMCHECK "Do cuda memcheck for testing?" OFF)
  if(FEAT_CUDAMEMCHECK)
    find_program(CUDA_MEMCHECK_EXE compute-sanitizer)
    set(CUDA_MEMCHECK_EXE ${CUDA_MEMCHECK_EXE} --leak-check full)
  endif(FEAT_CUDAMEMCHECK)

  option(FEAT_HAVE_HALFMATH "Support half precision floating point?" OFF)
  if(FEAT_HAVE_HALFMATH)
    target_include_directories(feat-dependency INTERFACE "${CUDAToolkit_INCLUDE_DIRS}")
  endif (FEAT_HAVE_HALFMATH)

endif(FEAT_HAVE_CUDA)

option(FEAT_HAVE_MPI "Enable MPI?" OFF)
if(FEAT_HAVE_MPI)
  if(FEAT_MPI_CXX)
    set(MPI_CXX_COMPILER ${FEAT_MPI_CXX})
    set(MPI_C_COMPILER ${FEAT_MPI_C})
  endif(FEAT_MPI_CXX)
  set(MPI_CXX_SKIP_MPICXX)
  find_package(MPI REQUIRED)
  execute_process(COMMAND "${MPIEXEC}" --version
    WORKING_DIRECTORY ${FEAT_SOURCE_DIR}
    OUTPUT_VARIABLE CMAKE_MPI_VERSION
    ERROR_VARIABLE CMAKE_MPI_VERSION
    )
  if(CMAKE_MPI_VERSION)
    string(REPLACE "\n" " " CMAKE_MPI_VERSION ${CMAKE_MPI_VERSION})
  endif(CMAKE_MPI_VERSION)

  ## TODO: also link against MPI_C ?
  target_link_libraries(feat-dependency INTERFACE MPI::MPI_CXX)

  #add_definitions( -DOMPI_SKIP_MPICXX )

  #set ctest runner to single process mpiexec execution
  set(FEAT_CTEST_RUNNER ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS})
endif(FEAT_HAVE_MPI)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  set(THREADS_PREFER_PTHREAD_FLAG TRUE)
  find_package(Threads REQUIRED)
  target_link_options(feat-dependency INTERFACE ${CMAKE_THREAD_LIBS_INIT})
endif()

option(FEAT_HAVE_OMP "FEAT OpenMP support available" OFF)
if(FEAT_HAVE_OMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(feat-dependency INTERFACE OpenMP::OpenMP_CXX)
endif(FEAT_HAVE_OMP)

option(FEAT_VALGRIND "Enable valgrind ctest support" OFF) #this should be done thorugh internal methods...
if(FEAT_VALGRIND)
  find_program(VALGRIND_EXE valgrind)
  set (VALGRIND_EXE ${VALGRIND_EXE} -v --suppressions=${FEAT_SOURCE_DIR}/valgrind.supp --fullpath-after=${FEAT_SOURCE_DIR} --leak-check=full --partial-loads-ok=yes)
  #--track-origins=yes
  #--keep-stacktraces=alloc-and-free
endif(FEAT_VALGRIND)

# our custom target to build (but not execute) all tests in all directories
add_custom_target( tests )

# custom target to build (but not execute) all cuda tests in all directories
if(FEAT_HAVE_CUDA)
  add_custom_target( cuda_tests )
endif()

# fetch contents outsourced to thirdparty subdirectory
add_subdirectory( thirdparty )
add_subdirectory( kernel )
add_subdirectory( test_system )
add_subdirectory( tools )
add_subdirectory( applications )
add_subdirectory( area51 )
add_subdirectory( benchmarks )
add_subdirectory( tutorials )
add_subdirectory( control )


if(Boost_FOUND)
message(STATUS "Boost Found ${Boost_FOUND}")
endif()
# store all variables describing the build in a custom header file
# to use later on in the real code
configure_file(${FEAT_SOURCE_DIR}/feat_config.hpp.in ${FEAT_BINARY_DIR}/feat_config.hpp)
