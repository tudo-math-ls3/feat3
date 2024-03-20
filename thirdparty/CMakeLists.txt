# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# add custom 'thirdparty' target
add_custom_target(thirdparty)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_ALGLIB)
  file(GLOB alglib-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ALGLIB}/cpp/src/*.cpp")

  add_library(thirdparty-alglib STATIC ${alglib-list})
  add_dependencies(thirdparty thirdparty-alglib)
  if (FEAT_COTIRE)
    cotire(thirdparty-alglib)
  endif (FEAT_COTIRE)
  target_compile_options(thirdparty-alglib PRIVATE -w)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(thirdparty-alglib PRIVATE -fno-lto)
    #set_property(TARGET thirdparty-alglib APPEND_STRING PROPERTY LINK_FLAGS " -fno-lto") #probably not need since no linking of a static library file...
  endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
endif (FEAT_HAVE_ALGLIB)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_FPARSER)
  file(GLOB fparser-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_FPARSER}/*.cc")

  # The 'fpoptimizer.cc' can not be compiled with clang++, so we have to disable
  # it by adding the corresponding pre-processor define
  if (${FEAT_COMPILER_ID} STREQUAL "clang")
    add_definitions(-DFP_NO_SUPPORT_OPTIMIZER)
  endif (${FEAT_COMPILER_ID} STREQUAL "clang")

  add_library(thirdparty-fparser STATIC ${fparser-list})
  add_dependencies(thirdparty thirdparty-fparser)
  target_compile_options(thirdparty-fparser PRIVATE -w)
endif (FEAT_HAVE_FPARSER)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_UMFPACK)
  file(GLOB suite_sparse-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/SuiteSparse_config/*.c")
    ##nice solution: swapp out the .c files to incorporate all variable assembly flags
  #always define our Source biuld through these wrappers
  # if(WIN32)
    file(GLOB amd-list "${FEAT_SOURCE_DIR}/thirdparty/SuiteSparseWrapperFiles/AMD/SourceWrappers/*.c")
    file(GLOB umfpack-list "${FEAT_SOURCE_DIR}/thirdparty/SuiteSparseWrapperFiles/UMFPACK/SourceWrappers/*.c")


  # else(WIN32)
  #   file(GLOB amd-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/AMD/Source/*.c")
  #   file(GLOB umfpack-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/UMFPACK/Source/*.c")
  # endif(WIN32)
  add_library(thirdparty-umfpack STATIC ${suite_sparse-list} ${amd-list} ${umfpack-list})
  add_dependencies(thirdparty thirdparty-umfpack)
  if (FEAT_COTIRE)
    cotire(thirdparty-umfpack)
  endif (FEAT_COTIRE)


  if (FEAT_HAVE_MKL AND FEAT_HAVE_OMP)
    target_compile_definitions(thirdparty-umfpack PRIVATE NCHOLMOD NTIMER LONGBLAS="long long")
    target_compile_options(thirdparty-umfpack PRIVATE -w)
    if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"))
      target_link_libraries(thirdparty-umfpack -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl)
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      target_link_libraries(thirdparty-umfpack -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl)
    else ()
      message(WARNING "compiler ${CMAKE_CXX_COMPILER_ID} not know in mkl library selection, using serial mkl library!")
      target_link_libraries(thirdparty-ufmpack -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl)
    endif ((${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"))
  elseif (FEAT_HAVE_MKL AND NOT FEAT_HAVE_OMP)
    target_compile_definitions(thirdparty-umfpack PRIVATE NCHOLMOD NTIMER LONGBLAS="long long")
    target_compile_options(thirdparty-umfpack PRIVATE -w)
    target_link_libraries(thirdparty-umfpack
       -Wl,--start-group libmkl_intel_ilp64.a libmkl_core.a libmkl_sequential.a -Wl,--end-group -lm -ldl)
    if (NOT FEAT_HAVE_MPI OR NOT "${MPI_CXX_LIBRARIES}" MATCHES ".*libpthread.*")
       target_link_libraries(thirdparty-umfpack -lpthread)
    endif (NOT FEAT_HAVE_MPI OR NOT "${MPI_CXX_LIBRARIES}" MATCHES ".*libpthread.*")
  else (FEAT_HAVE_MKL AND FEAT_HAVE_OMP)
    target_compile_definitions(thirdparty-umfpack PRIVATE NCHOLMOD NTIMER NBLAS)
    target_compile_options(thirdparty-umfpack PRIVATE -w)
  endif (FEAT_HAVE_MKL AND FEAT_HAVE_OMP)

  target_include_directories(thirdparty-umfpack PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/SuiteSparse_config/")
  target_include_directories(thirdparty-umfpack PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/AMD/Include/")
  target_include_directories(thirdparty-umfpack PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUITESPARSE}/UMFPACK/Include/")
endif (FEAT_HAVE_UMFPACK)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_PARMETIS)
  file(GLOB gklib-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/GKlib/*.c")
  file(GLOB metis-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/libmetis/*.c")
  add_library(thirdparty-metis STATIC ${gklib-list} ${metis-list})
  add_dependencies(thirdparty thirdparty-metis)
  if (FEAT_COTIRE)
    cotire(thirdparty-metis)
  endif (FEAT_COTIRE)
  target_compile_options(thirdparty-metis PRIVATE -w)
#  target_link_options(thirdparty-metis PRIVATE "LINKER:-w") #useless since static libraries do not use a linker
  target_include_directories(thirdparty-metis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/GKlib/")
  target_include_directories(thirdparty-metis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/include")
  target_include_directories(thirdparty-metis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/libmetis/")

  file(GLOB parmetis-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/libparmetis/*.c")
  add_library(thirdparty-parmetis STATIC ${parmetis-list})
  add_dependencies(thirdparty thirdparty-parmetis)
  if (FEAT_COTIRE)
    cotire(thirdparty-parmetis)
  endif (FEAT_COTIRE)
  target_compile_options(thirdparty-parmetis PRIVATE -w)
#  target_link_options(thirdparty-parmetis PRIVATE "LINKER:-w") #useless since static libraries do not use a linker
  target_include_directories(thirdparty-parmetis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/GKlib/")
  target_include_directories(thirdparty-parmetis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/metis/include/")
  target_include_directories(thirdparty-parmetis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/libparmetis")
  target_include_directories(thirdparty-parmetis PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_PARMETIS}/include")
  target_link_libraries(thirdparty-parmetis PRIVATE thirdparty-metis)
endif (FEAT_HAVE_PARMETIS)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_ZLIB)
  file(GLOB zlib-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZLIB}/*.c")
  add_library(thirdparty-zlib STATIC ${zlib-list})
  add_dependencies(thirdparty thirdparty-zlib)
  if (FEAT_COTIRE)
    cotire(thirdparty-zlib)
  endif (FEAT_COTIRE)
  target_compile_options(thirdparty-zlib PRIVATE -w)
#  target_link_options(thirdparty-zlib PRIVATE "LINKER:-w") #useless since static libraries do not use a linker
  #if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  #  set_property(TARGET thirdparty-zlib APPEND_STRING PROPERTY COMPILE_FLAGS " -fno-lto")
  #  set_property(TARGET thirdparty-zlib APPEND_STRING PROPERTY LINK_FLAGS " -fno-lto")
  #endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
endif (FEAT_HAVE_ZLIB)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_TRIANGLE)
  # only this one source file is required
  file(GLOB triangle-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_TRIANGLE}/triangle.c")

  add_library(thirdparty-triangle STATIC ${triangle-list})
  add_dependencies(thirdparty thirdparty-triangle)
  if (FEAT_COTIRE)
    cotire(thirdparty-triangle)
  endif (FEAT_COTIRE)
  target_compile_definitions(thirdparty-triangle PRIVATE ANSI_DECLARATORS TRILIBRARY)
  target_compile_options(thirdparty-triangle PRIVATE -w)
#  target_link_options(thirdparty-triangle PRIVATE "LINKER:-w")
endif (FEAT_HAVE_TRIANGLE)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_HYPRE)
  file(READ filelist_hypre hypre-list-in)
  string(REGEX REPLACE "\n" ";" hypre-list-in "${hypre-list-in}")
  foreach (hypre-list-item ${hypre-list-in})
    set(hypre-list "${hypre-list}" "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/${hypre-list-item}")
  endforeach ()
  add_library (thirdparty-hypre STATIC ${hypre-list})
  add_dependencies(thirdparty thirdparty-hypre)
  configure_file(${FEAT_SOURCE_DIR}/thirdparty/config_hypre.h.in ${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/HYPRE_config.h @ONLY)
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_ls")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/multivector")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/krylov")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/distributed_matrix")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/matrix_matrix")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/IJ_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_block_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/parcsr_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/seq_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/sstruct_ls")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/sstruct_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/struct_ls")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/struct_mv")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/utilities")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/blas")
  target_include_directories(thirdparty-hypre PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_HYPRE}/src/lapack")
  if (FEAT_HAVE_MPI)
    if (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
      target_compile_definitions(thirdparty-hypre PRIVATE HYPRE_HAVE_MPI)
      target_compile_definitions(thirdparty-hypre PRIVATE WIN32)
      target_compile_options(thirdparty-hypre PRIVATE -w)
    else (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
      target_compile_definitions(thirdparty-hypre PRIVATE HYPRE_HAVE_MPI)
      target_compile_options(thirdparty-hypre PRIVATE -w)
    endif (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
  else (FEAT_HAVE_MPI)
    if (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
      target_compile_options(thirdparty-hypre PRIVATE -w)
      target_compile_definitions(thirdparty-hypre PRIVATE HYPRE_SEQUENTIAL)
      target_compile_definitions(thirdparty-hypre PRIVATE WIN32)
    else (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
      target_compile_options(thirdparty-hypre PRIVATE -w)
      target_compile_definitions(thirdparty-hypre PRIVATE HYPRE_SEQUENTIAL)
    endif (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
  endif (FEAT_HAVE_MPI)
#  target_link_options(thirdparty-hypre PRIVATE "LINKER:-w")
endif (FEAT_HAVE_HYPRE)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_DEATH_HANDLER)
  add_library(thirdparty-deathhandler STATIC "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_DEATHHANDLER}/DeathHandler-master/death_handler.cc")
  add_dependencies(thirdparty thirdparty-deathhandler)
  target_compile_options(thirdparty-deathandler PRIVATE -w)
#  target_link_options(thirdparty-deathhandler PRIVATE "LINKER:-w")
endif (FEAT_HAVE_DEATH_HANDLER)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_ZFP)
  file(GLOB zfp-list "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZFP}/src/*.c")

  add_library(thirdparty-zfp STATIC ${zfp-list})
  add_dependencies(thirdparty thirdparty-zfp)
  if (FEAT_COTIRE)
    cotire(thirdparty-zfp)
  endif (FEAT_COTIRE)
  target_include_directories(thirdparty-zfp PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZFP}/src")
  target_include_directories(thirdparty-zfp PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZFP}/include")
  target_compile_options(thirdparty-zfp PRIVATE -w)
#  target_link_options(thirdparty-zfp PRIVATE "LINKER:-w")
endif (FEAT_HAVE_ZFP)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_SUPERLU_DIST)
  file(READ filelist_superlu superlu-list-in)
  string(REGEX REPLACE "\n" ";" superlu-list-in "${superlu-list-in}")
  foreach (superlu-list-item ${superlu-list-in})
    set(superlu-list "${superlu-list}" "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUPERLU}/${superlu-list-item}")
  endforeach ()
  add_library(thirdparty-superlu STATIC ${superlu-list})
  add_dependencies(thirdparty thirdparty-superlu)
  configure_file(${FEAT_SOURCE_DIR}/thirdparty/superlu_dist_config.h.in ${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUPERLU}/SRC/superlu_dist_config.h @ONLY)
  target_include_directories(thirdparty-superlu PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_SUPERLU}/SRC")
  set_property(TARGET thirdparty-superlu PROPERTY C_STANDARD 99)
  if (FEAT_HAVE_MPI)
    target_compile_options(thirdparty-superlu PRIVATE -w)
  endif (FEAT_HAVE_MPI)
#  target_link_options(thirdparty-superlu PRIVATE "LINKER:-w") #For link options, use LINKER generation to automatically transform to the correct linker call for compiler
endif (FEAT_HAVE_SUPERLU_DIST)

########################################################################################################################
########################################################################################################################
########################################################################################################################

if (FEAT_HAVE_ZOLTAN)
  file(READ filelist_zoltan zoltan-list-in)
  string(REGEX REPLACE "\n" ";" zoltan-list-in "${zoltan-list-in}")
  foreach (zoltan-list-item ${zoltan-list-in})
    set(zoltan-list "${zoltan-list}" "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/${zoltan-list-item}")
  endforeach ()
  add_library(thirdparty-zoltan STATIC ${zoltan-list})
  add_dependencies(thirdparty thirdparty-zoltan)
  configure_file(${FEAT_SOURCE_DIR}/thirdparty/config_zoltan.h.in ${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/include/Zoltan_config.h @ONLY)
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/include")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/all")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/Utilities/shared")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/Utilities/Timer")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/zz")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/coloring")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/graph")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/ha")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/hier")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/hsfc")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/lb")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/matrix")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/order")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/par")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/params")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/phg")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/rcb")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/reftree")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/simple")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/timer")
  target_include_directories(thirdparty-zoltan PRIVATE "${FEAT_SOURCE_DIR}/thirdparty/${FEAT_DIRNAME_ZOLTAN}/src/tpls")
  if (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
   #define __STDC__ since clang for some reason does not "confirm" to the standard,
   #but in this case this flag should be set to let zoltan compile
    target_compile_options(thirdparty-zoltan PRIVATE -w)
    target_compile_definitions(thirdparty-zoltan PRIVATE __STDC__=1)
  else (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")) #todo: Does this make sense? Only for WIN32?
    target_compile_options(thirdparty-zoltan PRIVATE -w)
  endif (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
#  target_link_options(thirdparty-zoltan PRIVATE "LINKER:-w")
endif (FEAT_HAVE_ZOLTAN)
