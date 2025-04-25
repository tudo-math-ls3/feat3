if(NOT TARGET Alglib::Alglib)
  FetchContent_GetProperties(alglib)

  file(GLOB alglib-srcs ${alglib_SOURCE_DIR}/src/*.cpp)

  add_library(feat-alglib-extern STATIC ${alglib-srcs})
  target_include_directories(feat-alglib-extern INTERFACE "${alglib_SOURCE_DIR}/src/")
  target_compile_options(feat-alglib-extern PRIVATE -w)

  add_library(Alglib::Alglib ALIAS feat-alglib-extern)
endif()
