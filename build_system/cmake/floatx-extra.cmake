if(NOT TARGET FloatX::FloatX)
  FetchContent_GetProperties(floatx)

  add_library(feat-floatx-extern INTERFACE)
  target_include_directories(feat-floatx-extern INTERFACE "${floatx_SOURCE_DIR}/src")
  add_library(FloatX::FloatX ALIAS feat-floatx-extern)
endif()
