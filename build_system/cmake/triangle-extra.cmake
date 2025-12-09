if(NOT TARGET triangle::triangle)
  FetchContent_GetProperties(triangle)

  add_library(feat-triangle-extern STATIC ${triangle_SOURCE_DIR}/triangle.c)
  target_compile_definitions(feat-triangle-extern PRIVATE ANSI_DECLARATORS TRILIBRARY)
  add_library(triangle::triangle ALIAS feat-triangle-extern)

  # disable warnings for this TPL
  if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    string(REGEX REPLACE "/W[1|2|3|4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    target_compile_options(feat-triangle-extern PRIVATE /W0)
  else()
    target_compile_options(feat-triangle-extern PRIVATE -w)
  endif()
endif()
