if(NOT TARGET fparser::fparser)
  FetchContent_GetProperties(fparser)

  file(GLOB fparser-srcs ${fparser_SOURCE_DIR}/*.cc)
  add_library(feat-fparser-extern STATIC ${fparser-srcs})
  target_include_directories(feat-fparser-extern INTERFACE ${fparser_SOURCE_DIR})
  target_compile_definitions(feat-fparser-extern PRIVATE FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS FP_USE_THREAD_SAFE_EVAL)

  # The 'fpoptimizer.cc' can not be compiled with clang++, so we have to disable
  # it by adding the corresponding pre-processor define
  if (${FEAT_COMPILER_ID} STREQUAL "clang")
    target_compile_definitions(feat-fparser-extern PRIVATE FP_NO_SUPPORT_OPTIMIZER)
  endif (${FEAT_COMPILER_ID} STREQUAL "clang")

  # disable warnings for this TPL
  if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    string(REGEX REPLACE "/W[1|2|3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE "/W[1|2|3|4]" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    target_compile_options(feat-fparser-extern PRIVATE /W0 /wd4244)
  else()
    target_compile_options(feat-fparser-extern PRIVATE -w)
  endif()

  add_library(fparser::fparser ALIAS feat-fparser-extern)
endif()
