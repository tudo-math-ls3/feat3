if(NOT TARGET zlib::zlib)
  FetchContent_GetProperties(zlib)

  file(GLOB zlib-sources "${zlib_SOURCE_DIR}/*.c")

  add_library(feat-zlib-extern STATIC ${zlib-sources})
  target_include_directories(feat-zlib-extern PUBLIC ${zlib_SOURCE_DIR})
  if(NOT WIN32)
    target_compile_options(feat-zlib-extern PRIVATE -w)
    target_compile_definitions(feat-zlib-extern PRIVATE HAVE_UNISTD_H)
  endif()
  add_library(zlib::zlib ALIAS feat-zlib-extern)
endif()
