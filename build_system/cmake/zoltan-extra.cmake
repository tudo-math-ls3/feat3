if(NOT TARGET Zoltan::Zoltan)
  FetchContent_GetProperties(zoltan)

  # Build file list
  file(READ filelist_zoltan zoltan-list-in)
  string(REGEX REPLACE "\n" ";" zoltan-list-in "${zoltan-list-in}")
  foreach (zoltan-list-item ${zoltan-list-in})
    set(zoltan-list "${zoltan-list}" "${zoltan_SOURCE_DIR}/${zoltan-list-item}")
  endforeach ()

  add_library(feat-zoltan-extern STATIC ${zoltan-list})
  configure_file(${FEAT_SOURCE_DIR}/thirdparty/config_zoltan.h.in ${zoltan_SOURCE_DIR}/src/include/Zoltan_config.h @ONLY)
  target_include_directories(feat-zoltan-extern PUBLIC  "${zoltan_SOURCE_DIR}/src/include")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/all")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/Utilities/shared")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/Utilities/Timer")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/zz")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/coloring")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/graph")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/ha")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/hier")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/hsfc")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/lb")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/matrix")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/order")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/par")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/params")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/phg")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/rcb")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/reftree")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/simple")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/timer")
  target_include_directories(feat-zoltan-extern PRIVATE "${zoltan_SOURCE_DIR}/src/tpls")

  if (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
   #define __STDC__ since clang for some reason does not "confirm" to the standard,
   #but in this case this flag should be set to let zoltan compile
    target_compile_options(feat-zoltan-extern PRIVATE -w)
    target_compile_definitions(feat-zoltan-extern PRIVATE __STDC__=1)
  else (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")) #todo: Does this make sense? Only for WIN32?
    target_compile_options(feat-zoltan-extern PRIVATE -w)
  endif (WIN32 AND (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    #target_link_options(thirdparty-zoltan PRIVATE "LINKER:-w")

  find_package(MPI REQUIRED)
  target_link_libraries(feat-zoltan-extern PRIVATE MPI::MPI_C)

  add_library(Zoltan::Zoltan ALIAS feat-zoltan-extern)
endif()
