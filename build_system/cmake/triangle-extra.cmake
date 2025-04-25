if(NOT TARGET triangle::triangle)
  FetchContent_GetProperties(triangle)

  add_library(feat-triangle-extern STATIC ${triangle_SOURCE_DIR}/triangle.c)
  target_compile_definitions(feat-triangle-extern PRIVATE ANSI_DECLARATORS TRILIBRARY)
  target_compile_options(feat-triangle-extern PRIVATE -w)
  add_library(triangle::triangle ALIAS feat-triangle-extern)
endif()
