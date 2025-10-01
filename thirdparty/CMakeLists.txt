# This file is responsible for finding/downloading third-party dependencies and
# including them in the FEAT3 build.

# By default enabled dependencies are first searched for among the libraries
# installed on the system.  If the dependency can not be found that way, it
# will automatically be downloaded.  Change the value of
# FETCHCONTENT_TRY_FIND_PACKAGE_MODE to "NEVER" to disable this behaviour.
# Downloaded dependencies are stored in a "_deps" folder in the build-tree.

# Third-party dependencies supported by FEAT3:
# +----------------+------------------------+----------------------------+
# |   Dependency   |          Flag          |           Target           |
# +----------------+------------------------+----------------------------+
# | Umfpack        | FEAT_HAVE_UMFPACK      | SuiteSparse::UMFPACK       |
# | Hypre          | FEAT_HAVE_HYPRE        | HYPRE::HYPRE               |
# | ParMETIS       | FEAT_HAVE_PARMETIS     | ParMETIS::ParMETIS         |
# | fparser        | FEAT_HAVE_FPARSER      | fparser::fparser           |
# | CGAL           | FEAT_HAVE_CGAL         | CGAL::CGAL                 |
# | Alglib         | FEAT_HAVE_ALGLIB       | Alglib::Alglib             |
# | deathhandler   | FEAT_HAVE_DEATHHANDLER | deathhandler::deathhandler |
# | floatx         | FEAT_HAVE_FLOATX       | FloatX::FloatX             |
# | ZFP            | FEAT_HAVE_ZFP          | ZFP::ZFP                   |
# | Zlib           | FEAT_HAVE_ZLIB         | zlib::zlib                 |
# | Triangle       | FEAT_HAVE_TRIANGLE     | triangle::triangle         |
# | Zoltan         | FEAT_HAVE_ZOLTAN       | Zoltan::Zoltan             |
# | SuperLU (dist) | FEAT_HAVE_SUPERLU_DIST | SuperLU::SuperLU           |
# | cuDSS          | FEAT_HAVE_CUDSS        | cuDSS::cuDSS               |
# | pmp            | FEAT_HAVE_PMP          | pmp::pmp                   |
# +----------------+------------------------+----------------------------+

# first of all, load fetchContent modul
include(FetchContent)
include(PrintPackageConfig)
include(FindPackageOverrideHelper)
include(FEATGetTPL)
include(CMakeDependentOption)

# Enable output during dependency download
set(FETCHCONTENT_QUIET OFF CACHE BOOL "fetch content verboseness")

option(FEAT_NO_EXTERNAL_DOWNLOAD "disable thirdparty downloads" OFF)

# Directory for cache of TPL archives.
# feat_get_tpl will search for archives in this directory and download missing archives to it.
# Unset (or set to an empty string) to disable caching behaviour.
set(FEAT_TPL_CACHE_DIRECTORY "${PROJECT_SOURCE_DIR}/thirdparty/cache" CACHE PATH "cache for downloaded TPL archives" )

# Only show the FEAT_PREFER_EXTERNAL_TPL option if FEAT_NO_EXTERNAL_DOWNLOAD is not true
# Otherwise force it to true
cmake_dependent_option(FEAT_PREFER_EXTERNAL_TPL
  "prefer installed libraries over downloaded ones" ON
  "NOT FEAT_NO_EXTERNAL_DOWNLOAD" ON)

# MAKE_AVAIL_LIST will contain all third-party dependencies that should be enabled
set(MAKE_AVAIL_LIST "")
mark_as_advanced(MAKE_AVAIL_LIST)

# Declare sources, i.e. git repositories / tarball urls, for third-party
# dependencies and set any applicable settings For sources we do not control,
# dependecies are pinned to specific versions either via git tags or md5sums.
# That way dependencies can not change without us noticing.

# UPDATE_DISCONNECTED ON is set for all thirdparty dependencies. This means
# cmake will not try to update dependencies unless any details of the download
# or update steps have changed.  This is required for reliable patching of
# dependencies. Otherwise CMake would try to pull the newest changes from the
# git repositories, which would cause the configuring to fail if any of the
# PATCH_COMMANDS cause a git conflict. Because all our dependcies are pinned
# anyway, this isn't a loss of functionality.

# Add "-w" to the C/CXX flags to disable warnings for third-party dependencies.
# Save the current flags beforehand to be able to restore them after handling
# third-party dependencies.
set(CMAKE_CXX_FLAGS_BACKUP ${CMAKE_CXX_FLAGS})
set(CMAKE_C_FLAGS_BACKUP ${CMAKE_C_FLAGS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")

if(FEAT_TPL_CACHE_DIRECTORY)
  file(MAKE_DIRECTORY ${FEAT_TPL_CACHE_DIRECTORY})
endif()

option(FEAT_HAVE_UMFPACK "Use umfpack solver of suitesparse?" OFF)
if(FEAT_HAVE_UMFPACK)
  feat_get_tpl(
    PACKAGE_NAME UMFPACK
    VERSION 6.3.5
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.10.2.zip
    URL_HASH MD5=b54b521b41f1468e37e0c0199a85e598
    CONFIG
  )

  if(NOT UMFPACK_FOUND)
    ## set cmake variables for suitesparse
    # this should build everything necessary for umfpack
    set(SUITESPARSE_ENABLE_PROJECTS "umfpack" CACHE STRING "which projects to build?" FORCE)
    # build static libraries, are these options local?
    set(BUILD_STATIC_LIBS "ON" CACHE BOOL "build static libs" FORCE)
    # do not build shared libries
    set(BUILD_SHARED_LIBS "OFF" CACHE BOOL "build shared libs" FORCE)
    # i think we should disable cuda
    set(SUITESPARSE_USE_CUDA "OFF" CACHE BOOL "build with cuda support" FORCE)
    # do not use fortran compiler to determine blas function names
    # warning: this requires SUITESPARSE_C_TO_FORTRAN to be set correctly
    set(SUITESPARSE_USE_FORTRAN "OFF" CACHE BOOL "build with cuda support" FORCE)
    # Enable Suitesparse OpenMP only if FEAT OpenMP is enabled
    set(SUITESPARSE_USE_OPENMP ${FEAT_HAVE_OMP} CACHE BOOL "build with OpenMP support" FORCE)

    # TODO: Check if we can find out if the available blas is 32 or 64 bit
    set(SUITESPARSE_USE_64BIT_BLAS "ON" CACHE BOOL "use 64 bit blas" FORCE)

    # Disable CHOLDMOD tests
    set(BUILD_TESTING "OFF")
  endif()
endif()

option(FEAT_HAVE_HYPRE "Use hypre solvers?" OFF)
if(FEAT_HAVE_HYPRE)
  feat_get_tpl(
    PACKAGE_NAME HYPRE
    VERSION 2.33
    URL https://github.com/hypre-space/hypre/archive/refs/tags/v2.33.0.zip
    URL_HASH MD5=6fe2c8ef054836d88bb8fad9794c4db5
    PATCH_COMMAND_WINDOWS python ${CMAKE_CURRENT_LIST_DIR}/patches/hypre_win.py .
    SOURCE_SUBDIR src
    CONFIG
  )

  if(NOT HYPRE_FOUND)
    # Enable hypre MPI support if FEAT is build with MPI support
    set(HYPRE_ENABLE_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)
    set(HYPRE_HAVE_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)
    set(HYPRE_WITH_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)

    # Enable Hypre OpenMP support if FEAT is build with OpenMP support
    set(HYPRE_WITH_OPENMP ${FEAT_HAVE_OMP} CACHE BOOL "Build Hypre with OMP support" FORCE)
  endif()
endif()

option(FEAT_HAVE_PARMETIS "Enable ParMETIS" OFF)
if(FEAT_HAVE_PARMETIS)
  feat_get_tpl(
    PACKAGE_NAME parmetis
    VERSION 4.0.3
    URL https://github.com/otkmesser/ParMETIS/archive/refs/tags/v4.0.3.zip
    URL_HASH MD5=9bb15981d34adc3fd9e038bc7033388d
    PATCH_COMMAND_WINDOWS python ${CMAKE_CURRENT_LIST_DIR}/patches/parmetis_win.py .
    CONFIG
  )

endif(FEAT_HAVE_PARMETIS)

option(FEAT_HAVE_FPARSER "Enable FParser" OFF)
if(FEAT_HAVE_FPARSER)
  feat_get_tpl(
    PACKAGE_NAME fparser
    VERSION 4.5.2
    URL http://warp.povusers.org/FunctionParser/fparser4.5.2.zip
    URL_HASH MD5=ce805d2f827f25e61f266e176b8d338d
  )
endif()

option(FEAT_HAVE_BOOST "Enable Boost" OFF)
if(FEAT_HAVE_BOOST)
  feat_get_tpl(
    PACKAGE_NAME Boost
    VERSION 1.88
    URL https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.zip
    URL_HASH MD5=419f6a9273cb90d4f3f65ca0ae02cd00
  )

  if(NOT Boost_FOUND)
    # Only build libraries required for CGAL
    set(BOOST_INCLUDE_LIBRARIES graph heap logic)
    # Make boost create Boost::boost target
    set(BOOST_ENABLE_COMPATIBILITY_TARGETS ON)
  endif()
endif()

option(FEAT_HAVE_CGAL "Enable CGAL" OFF)
if(FEAT_HAVE_CGAL)
  feat_get_tpl(
    PACKAGE_NAME CGAL
    VERSION 6.0.1
    URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.zip
    URL_HASH MD5=f30c69bf83e294c44e936fa4ad2c2b2c
    SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
    CONFIG
  )
endif()

option(FEAT_HAVE_ALGLIB "Enable Alglib" OFF)
if(FEAT_HAVE_ALGLIB)
  feat_get_tpl(
    PACKAGE_NAME Alglib
    VERSION 4.04.0
    URL http://www.alglib.net/translator/re/alglib-4.04.0.cpp.gpl.zip
    URL_HASH MD5=c9a6b38c9b41ef1121fb0b33c8228394
  )
endif()

option(FEAT_HAVE_DEATH_HANDLER "Enable deathhandler" OFF)
if(FEAT_HAVE_DEATH_HANDLER)
  feat_get_tpl(
    PACKAGE_NAME death_handler
    VERSION 1.0
    URL https://github.com/vmarkovtsev/DeathHandler/archive/refs/heads/master.zip
    URL_HASH MD5=c1e443dbd6b1e5a6e0cb072a5726e9a0
  )
endif()

option(FEAT_HAVE_FLOATX "Enable FloatX" OFF)
if(FEAT_HAVE_FLOATX)
  feat_get_tpl(
    PACKAGE_NAME FloatX
    VERSION 1.0
    URL https://github.com/oprecomp/FloatX/archive/refs/heads/master.zip
    URL_HASH MD5=e8bd646b92de39ee2404ecafd08941de
    SOURCE_SUBDIR Non-existing
  )
endif()

option(FEAT_HAVE_ZFP "Enable zfp" OFF)
if(FEAT_HAVE_ZFP)
  feat_get_tpl(
    PACKAGE_NAME zfp
    VERSION 1.0.0
    URL https://github.com/LLNL/zfp/releases/download/1.0.0/zfp-1.0.0.zip
    URL_HASH MD5=28c1a0a53931973b18eb4b3d3b7c5f63
    CONFIG
  )
  if(NOT zfp_FOUND)
    set(ZFP_WITH_OPENMP ${FEAT_HAVE_OMP} CACHE BOOL "Enable OpenMP for ZFP" FORCE)
    # Disable zfp tests
    set(BUILD_TESTING "OFF")
 endif()
endif()

option(FEAT_HAVE_ZLIB "Enable zlib" OFF)
if(FEAT_HAVE_ZLIB)
  feat_get_tpl(
    PACKAGE_NAME zlib
    VERSION 1.2.13
    URL https://www.zlib.net/fossils/zlib-1.2.13.tar.gz
    URL_HASH MD5=9b8aa094c4e5765dabf4da391f00d15c
    SOURCE_SUBDIR Non-existing
  )
endif()

option(FEAT_HAVE_TRIANGLE "Enable triangle" OFF)
if(FEAT_HAVE_TRIANGLE)
  feat_get_tpl(
    PACKAGE_NAME triangle
    VERSION 1.0
    URL "http://www.netlib.org/voronoi/triangle.zip"
    URL_HASH MD5=10aff8d7950f5e0e2fb6dd2e340be2c9
  )
endif()

option(FEAT_HAVE_ZOLTAN "Enable zoltan" OFF)
if(FEAT_HAVE_ZOLTAN)
  feat_get_tpl(
    PACKAGE_NAME zoltan
    VERSION 3.901
    URL "https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz"
    URL_HASH MD5=a951ef9b1529fed3a325d57e0be21734
    PATCH_COMMAND_WINDOWS python ${CMAKE_CURRENT_LIST_DIR}/patches/zoltan_win.py .
    SOURCE_SUBDIR Non-Existing
  )
endif()

option(FEAT_HAVE_SUPERLU_DIST "Enable SuperLU (dist)" OFF)
if(FEAT_HAVE_SUPERLU_DIST)
  feat_get_tpl(
    PACKAGE_NAME superlu
    VERSION 9.1.0
    URL "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v9.1.0.tar.gz"
    URL_HASH MD5=f8559a94ee64c8b70ebd79b65576d08d
    PATCH_COMMAND_WINDOWS python ${CMAKE_CURRENT_LIST_DIR}/patches/superlu_win.py .
    PATCH_COMMAND_LINUX python ${CMAKE_CURRENT_LIST_DIR}/patches/superlu.py .
    SOURCE_SUBDIR SRC
    EXCLUDE_FROM_ALL
  )

  if(NOT superlu_FOUND)
    set(VERSION_MAJOR 1)
    set(SUPERLU_DIST_ENABLE_DOUBLE ON)
    set(SUPERLU_DIST_ENABLE_OMP ${FEAT_HAVE_OMP})

    set(SUPERLU_DIST_TPL_ENABLE_INTERNAL_BLASLIB ON)
    set(SUPERLU_DIST_TPL_ENABLE_PARMETIS ${FEAT_HAVE_PARMETIS})
    set(SUPERLU_DIST_ENABLE_COMPLEX_16 OFF)
    set(SUPERLU_DIST_ENABLE_SINGLE OFF)

    # appearently, we need to set the index size to 64 if combined with parmetis...
    # todo: get this from parmetis...
    if(SUPERLU_DIST_TPL_ENABLE_PARMETIS)
      set(XSDK_INDEX_SIZE 64)
    endif()
  endif()
endif()

option(FEAT_HAVE_CUDSS "Enable cuDSS" OFF)
if(FEAT_HAVE_CUDSS)
  if(WIN32)
    set(cudss_url "https://developer.download.nvidia.com/compute/cudss/redist/libcudss/windows-x86_64/libcudss-windows-x86_64-0.5.0.16_cuda12-archive.zip")
    set(cudss_hash "2ea2f64901db68c0e0e7f419ab60bef1")
  else()
    set(cudss_url "https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.5.0.16_cuda12-archive.tar.xz")
    set(cudss_hash "cb438f906120c01d062bd5db40244cd2")
  endif()

  feat_get_tpl(
    PACKAGE_NAME cudss
    VERSION 0.5.0.16
    URL ${cudss_url}
    URL_HASH MD5=${cudss_hash}
  )
endif()

option(FEAT_HAVE_PMP "Enable pmp" OFF)
if(FEAT_HAVE_PMP)
  feat_get_tpl(
    PACKAGE_NAME pmp
    VERSION 3.0.0
    URL https://github.com/pmp-library/pmp-library/archive/refs/tags/3.0.0.zip
    URL_HASH MD5=7b7f9ce07a7a687c9d78a6583cf64a2c
    SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
    CONFIG
  )
endif()

# after all target are declared, call now make available
if(MAKE_AVAIL_LIST)
  message(STATUS "Fetching the following third-party dependencies:  ${MAKE_AVAIL_LIST}")
  foreach(tpl IN LISTS MAKE_AVAIL_LIST)
    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...")
    message(STATUS "-------------------------------------------------------------------")
    FetchContent_MakeAvailable(${tpl})

    message(STATUS "-------------------------------------------------------------------")
    message(STATUS "- Configuring ${tpl}...done")
    message(STATUS "-------------------------------------------------------------------")
  endforeach()
endif()

# Restore previous C/CXX flags to re-enable warnings
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_BACKUP})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_BACKUP})
