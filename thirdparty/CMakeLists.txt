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

# Should dependecies be looked for among the system libraries?
# set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE OPT_IN CACHE STRING "fetch content mode")
# Enable output during dependency download
set(FETCHCONTENT_QUIET OFF CACHE BOOL "fetch content verboseness")
set(FEAT_NO_EXTERNAL_DOWNLOAD OFF CACHE BOOL "no thirdparty download allowed")

set(FEAT_PREFER_EXTERNAL_TPL ON CACHE BOOL "try to find external libraries before fetching?")
if(FEAT_NO_EXTERNAL_DOWNLOAD)
  message(STATUS "External tpl implicitly required")
  set(FEAT_PREFER_EXTERNAL_TPL ON CACHE BOOL "try to find external libraries before fetching?" FORCE)
endif()


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

option(FEAT_HAVE_UMFPACK "Use umfpack solver of suitesparse?" OFF)
if(FEAT_HAVE_UMFPACK)
  # TODO: instead use value of FETCHCONTENT_TRY_FIND_PACKAGE_MODE?
  if(FEAT_PREFER_EXTERNAL_TPL OR UMFPACK_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(UMFPACK CONFIG REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(UMFPACK CONFIG GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT UMFPACK_FOUND)
    FetchContent_Declare(
      # suitesparse
      UMFPACK
      GIT_REPOSITORY      https://github.com/DrTimothyAldenDavis/SuiteSparse
      GIT_TAG             v7.10.2
      GIT_PROGRESS ON
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS NAMES UMFPACK GLOBAL
      OVERRIDE_FIND_PACKAGE
    )

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

    # list(APPEND MAKE_AVAIL_LIST suitesparse)
    list(APPEND MAKE_AVAIL_LIST UMFPACK)
    find_package_override_helper(umfpack 7.10.2 AnyNewerVersion)
  else()
    print_package_info(UMFPACK)
  endif(NOT UMFPACK_FOUND)
endif(FEAT_HAVE_UMFPACK)

option(FEAT_HAVE_HYPRE "Use hypre solvers?" OFF)
if(FEAT_HAVE_HYPRE)
  if(FEAT_PREFER_EXTERNAL_TPL OR HYPRE_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(HYPRE CONFIG REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(HYPRE CONFIG GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT HYPRE_FOUND)
    if(WIN32)
      FetchContent_Declare(
        HYPRE
        GIT_REPOSITORY      https://github.com/hypre-space/hypre
        GIT_TAG             d3243638e7dee5208cb4ea643610fa9bdb74123b # v2.33.0
        GIT_PROGRESS        ON
        PATCH_COMMAND       python ${CMAKE_CURRENT_LIST_DIR}/patches/hypre_win.py .
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR       src # HYPRE's CMakeLists.txt is in roots src subdirectory, instead of the project root
        # FIND_PACKAGE_ARGS NAMES HYPRE GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    else()
      FetchContent_Declare(
        HYPRE
        GIT_REPOSITORY      https://github.com/hypre-space/hypre
        GIT_TAG             d3243638e7dee5208cb4ea643610fa9bdb74123b # v2.33.0
        GIT_PROGRESS        ON
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR       src # HYPRE's CMakeLists.txt is in roots src subdirectory, instead of the project root
        # FIND_PACKAGE_ARGS NAMES HYPRE GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    endif()

    list(APPEND MAKE_AVAIL_LIST HYPRE)
    find_package_override_helper(hypre 2.33 AnyNewerVersion)

    # Enable hypre MPI support if FEAT is build with MPI support
    set(HYPRE_ENABLE_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)
    set(HYPRE_HAVE_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)
    set(HYPRE_WITH_MPI ${FEAT_HAVE_MPI} CACHE BOOL "Build Hypre with MPI support" FORCE)

    # Enable Hypre OpenMP support if FEAT is build with OpenMP support
    set(HYPRE_WITH_OPENMP ${FEAT_HAVE_OMP} CACHE BOOL "Build Hypre with OMP support" FORCE)
  else()
    print_package_info(HYPRE)
  endif(NOT HYPRE_FOUND)
endif(FEAT_HAVE_HYPRE)

option(FEAT_HAVE_PARMETIS, "Enable ParMETIS" OFF)
if(FEAT_HAVE_PARMETIS)

  if(FEAT_PREFER_EXTERNAL_TPL OR parmetis_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(parmetis REQUIRED GLOBAL CONFIG)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(parmetis GLOBAL CONFIG)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT parmetis_FOUND)
    if(WIN32)
      FetchContent_Declare(
        parmetis
        GIT_REPOSITORY https://github.com/otkmesser/ParMETIS
        GIT_TAG main
        GIT_PROGRESS ON
        PATCH_COMMAND python ${CMAKE_CURRENT_LIST_DIR}/patches/parmetis_win.py .
        UPDATE_DISCONNECTED ON
        # FIND_PACKAGE_ARGS NAMES parmetis GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    else()
      FetchContent_Declare(
        parmetis
        GIT_REPOSITORY https://github.com/otkmesser/ParMETIS
        GIT_TAG main
        GIT_PROGRESS ON
        UPDATE_DISCONNECTED ON
        # FIND_PACKAGE_ARGS NAMES parmetis GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    endif()
    list(APPEND MAKE_AVAIL_LIST parmetis)
    find_package_override_helper(parmetis 4.0.3 AnyNewerVersion)
  else()
    print_package_info(parmetis)
  endif(NOT parmetis_FOUND)
endif(FEAT_HAVE_PARMETIS)

option(FEAT_HAVE_FPARSER "Enable FParser" OFF)
if(FEAT_HAVE_FPARSER)
  if(FEAT_PREFER_EXTERNAL_TPL OR fparser_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(fparser REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(fparser GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT fparser)
    FetchContent_declare(
      fparser
      URL http://warp.povusers.org/FunctionParser/fparser4.5.2.zip
      URL_HASH MD5=ce805d2f827f25e61f266e176b8d338d
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST fparser)
    find_package_override_helper(fparser 4.5.2 AnyNewerVersion)
  else()
    print_package_info(fparser)
  endif(NOT fparser)
endif(FEAT_HAVE_FPARSER)

option(FEAT_HAVE_BOOST "Enable Boost" OFF)
if(FEAT_HAVE_BOOST)

  if(FEAT_PREFER_EXTERNAL_TPL OR Boost_DIR)
    # TODO: Raise minimum boost to 1.82? See https://github.com/boostorg/cmake/tree/7f5336b3bf8067bb40da4e8b9940c133271e938a
    # before targets are not correctly set
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(Boost REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(Boost GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT Boost_FOUND)
    set(BOOST_INCLUDE_LIBRARIES graph heap logic)
    # Fetch Boost from the github repository. For some reason the github releases
    # contain a root CMakeLists.txt and the website releases do not.
    FetchContent_declare(
      Boost
      URL https://github.com/boostorg/boost/releases/download/boost-1.81.0/boost-1.81.0.tar.xz
      URL_HASH MD5=6cf0cdd797bca685910d527ae3c08cb3
      PATCH_COMMAND python ${CMAKE_CURRENT_LIST_DIR}/patches/boost.py .
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS CONFIG GLOBAL
      OVERRIDE_FIND_PACKAGE
      # EXCLUDE_FROM_ALL #TODO: probably a good idea
    )
    list(APPEND MAKE_AVAIL_LIST Boost)
    find_package_override_helper(boost 1.81 AnyNewerVersion)
  else()
    print_package_info(Boost)
  endif(NOT Boost_FOUND)
endif(FEAT_HAVE_BOOST)

option(FEAT_HAVE_CGAL "Enable CGAL" OFF)
if(FEAT_HAVE_CGAL)

  if(FEAT_PREFER_EXTERNAL_TPL OR CGAL_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(CGAL CONFIG REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(CGAL CONFIG GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT CGAL_FOUND)
    FetchContent_declare(
      CGAL
      URL https://github.com/CGAL/cgal/releases/download/v5.5.5/CGAL-5.5.5.tar.xz
      URL_HASH MD5=f17ec1a7667a6ebb28137dedcd327b74
      SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS NAMES CGAL GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST CGAL)
    find_package_override_helper(cgal 5.5.5 AnyNewerVersion)
  else()
    print_package_info(CGAL)
  endif(NOT CGAL_FOUND)
endif(FEAT_HAVE_CGAL)

option(FEAT_HAVE_ALGLIB "Enable Alglib" OFF)
if(FEAT_HAVE_ALGLIB)
  if(FEAT_PREFER_EXTERNAL_TPL OR Alglib_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(Alglib REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(Alglib GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT Alglib_FOUND)
    FetchContent_declare(
      Alglib
      URL http://www.alglib.net/translator/re/alglib-3.16.0.cpp.gpl.zip
      URL_HASH MD5=b12c3babddaa74ad634116a2bbed06a3
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS NAMES Alglib GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST Alglib)
    find_package_override_helper(alglib 3.16.0 AnyNewerVersion)
  else()
    print_package_info(Alglib)
  endif(NOT Alglib_FOUND)
endif(FEAT_HAVE_ALGLIB)

option(FEAT_HAVE_DEATH_HANDLER "Enable deathhandler" OFF)
if(FEAT_HAVE_DEATH_HANDLER)
  if(FEAT_PREFER_EXTERNAL_TPL OR death_handler_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(death_handler REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(death_handler GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT death_handler_FOUND)
    FetchContent_declare(
      death_handler
      GIT_REPOSITORY https://github.com/vmarkovtsev/DeathHandler
      GIT_TAG 911e4df429a89309b10111d76d232d4a7c2b5485
      UPDATE_DISCONNECTED ON
      GIT_PROGRESS ON
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST death_handler)
    find_package_override_helper(death_handler 1.0 AnyNewerVersion)
  else()
    print_package_info(death_handler)
  endif(NOT death_handler_FOUND)
endif()

option(FEAT_HAVE_FLOATX "Enable FloatX" OFF)
if(FEAT_HAVE_FLOATX)
  if(FEAT_PREFER_EXTERNAL_TPL OR FloatX_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(FloatX REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(FloatX GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT FloatX_FOUND)
    FetchContent_declare(
      FloatX
      GIT_REPOSITORY https://github.com/oprecomp/FloatX
      GIT_TAG 5ae867b16310864a6d7654a0253b0064289cf822
      GIT_PROGRESS ON
      UPDATE_DISCONNECTED ON
      SOURCE_SUBDIR Non-Existing
      # FIND_PACKAGE_ARGS NAMES FloatX GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST FloatX)
    find_package_override_helper(FloatX 1.0 AnyNewerVersion)
  else()
    print_package_info(FloatX)
  endif(NOT FloatX_FOUND)
endif()

option(FEAT_HAVE_ZFP "Enable zfp" OFF)
if(FEAT_HAVE_ZFP)
  if(FEAT_PREFER_EXTERNAL_TPL OR zfp_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zfp CONFIG REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zfp CONFIG GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT zfp_FOUND)
    FetchContent_declare(
      zfp
      GIT_REPOSITORY https://github.com/LLNL/zfp
      GIT_TAG f39af72648a2aeb88e9b2cca8c64f51b493ad5f4 #v1.0.0
      GIT_PROGRESS ON
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )

    set(ZFP_WITH_OPENMP ${FEAT_HAVE_OMP} CACHE BOOL "Enable OpenMP for ZFP" FORCE)

    # Disable zfp tests
    set(BUILD_TESTING "OFF")

    list(APPEND MAKE_AVAIL_LIST zfp)
    find_package_override_helper(zfp 1.0.0 AnyNewerVersion)
  else()
    print_package_info(zfp)
  endif(NOT zfp_FOUND)
endif()

option(FEAT_HAVE_ZLIB "Enable zlib" OFF)
if(FEAT_HAVE_ZLIB)
  if(FEAT_PREFER_EXTERNAL_TPL OR zlib_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zlib REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zlib GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT zlib_FOUND)
    FetchContent_declare(
      zlib
      URL https://www.zlib.net/fossils/zlib-1.2.13.tar.gz
      URL_HASH MD5=9b8aa094c4e5765dabf4da391f00d15c
      UPDATE_DISCONNECTED ON
      # zlibs default CMakeLists.txt adds example binaries to the build. Ignore it.
      SOURCE_SUBDIR Non-Existing
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST zlib)
    find_package_override_helper(zlib 1.2.13 AnyNewerVersion)
  else()
    print_package_info(zlib)
  endif(NOT zlib_FOUND)
endif()

option(FEAT_HAVE_TRIANGLE "Enable triangle" OFF)
if(FEAT_HAVE_TRIANGLE)
  if(FEAT_PREFER_EXTERNAL_TPL OR triangle_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(triangle REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(triangle GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT triangle_FOUND)
    FetchContent_declare(
      triangle
      URL "http://www.netlib.org/voronoi/triangle.zip"
      URL_HASH MD5=10aff8d7950f5e0e2fb6dd2e340be2c9
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )
    list(APPEND MAKE_AVAIL_LIST triangle)
    find_package_override_helper(triangle 1.0 AnyNewerVersion)
  else()
    print_package_info(triangle)
  endif(NOT triangle_FOUND)
endif()

option(FEAT_HAVE_ZOLTAN "Enable zoltan" OFF)
if(FEAT_HAVE_ZOLTAN)
  if(FEAT_PREFER_EXTERNAL_TPL OR zoltan_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zoltan REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(zoltan GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT zoltan_FOUND)
    if(WIN32)
      FetchContent_declare(
        zoltan
        URL "https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz"
        URL_HASH MD5=a951ef9b1529fed3a325d57e0be21734
        PATCH_COMMAND python ${CMAKE_CURRENT_LIST_DIR}/patches/zoltan_win.py .
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR Non-Existing
        # FIND_PACKAGE_ARGS GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    else()
      FetchContent_declare(
        zoltan
        URL "https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz"
        URL_HASH MD5=a951ef9b1529fed3a325d57e0be21734
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR Non-Existing
        # FIND_PACKAGE_ARGS GLOBAL
        OVERRIDE_FIND_PACKAGE
      )
    endif()
    list(APPEND MAKE_AVAIL_LIST zoltan)
    find_package_override_helper(zoltan 3.901 AnyNewerVersion)
  else()
    print_package_info(zoltan)
  endif(NOT zoltan_FOUND)
endif()

option(FEAT_HAVE_SUPERLU_DIST "Enable SuperLU (dist)" OFF)
if(FEAT_HAVE_SUPERLU_DIST)
  if(FEAT_PREFER_EXTERNAL_TPL OR superlu_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(superlu REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(superlu GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT superlu_FOUND)
    if(WIN32)
      ## probably broken right now, since we need to patch and use the new cmakelists...
      FetchContent_declare(
        superlu
        URL "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v9.1.0.tar.gz"
        URL_HASH MD5=f8559a94ee64c8b70ebd79b65576d08d
        PATCH_COMMAND python ${CMAKE_CURRENT_LIST_DIR}/patches/superlu_win.py .
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR SRC
        OVERRIDE_FIND_PACKAGE
        EXCLUDE_FROM_ALL
      )
    else()
      FetchContent_declare(
        superlu
        URL "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v9.1.0.tar.gz"
        URL_HASH MD5=f8559a94ee64c8b70ebd79b65576d08d
        PATCH_COMMAND python ${CMAKE_CURRENT_LIST_DIR}/patches/superlu.py .
        UPDATE_DISCONNECTED ON
        SOURCE_SUBDIR SRC
        OVERRIDE_FIND_PACKAGE
        EXCLUDE_FROM_ALL
      )
    endif()

    set(VERSION_MAJOR 1)
    set(SUPERLU_DIST_ENABLE_DOUBLE ON)
    set(SUPERLU_DIST_ENABLE_OMP ${FEAT_HAVE_OMP})

    set(SUPERLU_DIST_TPL_ENABLE_INTERNAL_BLASLIB ON)
    set(SUPERLU_DIST_TPL_ENABLE_PARMETISLIB ${FEAT_HAVE_PARMETIS})
    set(SUPERLU_DIST_ENABLE_COMPLEX_16 OFF)
    set(SUPERLU_DIST_ENABLE_SINGLE OFF)
    # appearently, we need to set the index size to 64 if combined with parmetis...
    # todo: get this from parmetis...
    if(SUPERLU_DIST_TPL_ENABLE_PARMETISLIB)
      set(XSDK_INDEX_SIZE 64)
    endif()

    list(APPEND MAKE_AVAIL_LIST superlu)
    find_package_override_helper(superlu 9.1.0 AnyNewerVersion)
  else()
    print_package_info(superlu)
  endif(NOT superlu_FOUND)
endif()

option(FEAT_HAVE_CUDSS "Enable cuDSS" OFF)
if(FEAT_HAVE_CUDSS)
  if(FEAT_PREFER_EXTERNAL_TPL OR cudss_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(cudss REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(cudss GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT cudss_FOUND)
    if(WIN32)
      set(cudss_url "https://developer.download.nvidia.com/compute/cudss/redist/libcudss/windows-x86_64/libcudss-windows-x86_64-0.5.0.16_cuda12-archive.zip")
      set(cudss_hash "2ea2f64901db68c0e0e7f419ab60bef1")
    else()
      set(cudss_url "https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.5.0.16_cuda12-archive.tar.xz")
      set(cudss_hash "cb438f906120c01d062bd5db40244cd2")
    endif()

    FetchContent_declare(
      cudss
      URL ${cudss_url}
      URL_HASH MD5=${cudss_hash}
      UPDATE_DISCONNECTED ON
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )

    list(APPEND MAKE_AVAIL_LIST cudss)
    find_package_override_helper(cudss 0.5.0.16 AnyNewerVersion)
  else()
    print_package_info(cudss)
  endif(NOT cudss_FOUND)
endif()

option(FEAT_HAVE_PMP "Enable pmp" OFF)
if(FEAT_HAVE_PMP)
  if(FEAT_PREFER_EXTERNAL_TPL OR pmp_DIR)
    if(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(pmp CONFIG REQUIRED GLOBAL)
    else(FEAT_NO_EXTERNAL_DOWNLOAD)
      find_package(pmp CONFIG GLOBAL)
    endif(FEAT_NO_EXTERNAL_DOWNLOAD)
  endif()

  if(NOT pmp_FOUND)
    FetchContent_declare(
      pmp
      GIT_REPOSITORY https://github.com/pmp-library/pmp-library
      GIT_TAG 3.0.0
      GIT_PROGRESS ON
      UPDATE_DISCONNECTED ON
      SOURCE_SUBDIR Non-Existing # Use non-existing source dir to disable add_subdirectory call of MakeAvailable
      # FIND_PACKAGE_ARGS GLOBAL
      OVERRIDE_FIND_PACKAGE
    )

    list(APPEND MAKE_AVAIL_LIST pmp)
    find_package_override_helper(pmp 3.0.0 AnyNewerVersion)
  else()
    print_package_info(pmp)
  endif(NOT pmp_FOUND)
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
