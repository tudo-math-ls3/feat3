if(NOT TARGET zlib::zlib)
  FetchContent_GetProperties(zlib)

  file(GLOB zlib-sources "${zlib_SOURCE_DIR}/*.c")
  add_library(feat-zlib-extern STATIC ${zlib-sources})
  target_compile_options(feat-zlib-extern PRIVATE -w)
  target_compile_definitions(feat-zlib-extern PRIVATE HAVE_UNISTD_H)
  add_library(zlib::zlib ALIAS feat-zlib-extern)
endif()
