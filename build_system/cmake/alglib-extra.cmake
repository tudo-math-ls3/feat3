if(NOT TARGET Alglib::Alglib)
  FetchContent_GetProperties(alglib)

  file(GLOB alglib-srcs ${alglib_SOURCE_DIR}/src/*.cpp)

  add_library(feat-alglib-extern STATIC ${alglib-srcs})
  target_include_directories(feat-alglib-extern INTERFACE "${alglib_SOURCE_DIR}/src/")

  # disable warnings for this TPL
  if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    string(REGEX REPLACE "/W[1|2|3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE "/W[1|2|3|4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    target_compile_options(feat-alglib-extern PRIVATE /W0)
  else()
    target_compile_options(feat-alglib-extern PRIVATE -w)
  endif()

  add_library(Alglib::Alglib ALIAS feat-alglib-extern)
endif()
