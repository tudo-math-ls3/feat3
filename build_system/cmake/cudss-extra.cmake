if(NOT TARGET cuDSS::cuDSS)
  FetchContent_GetProperties(cudss)

  add_library(feat-cudss-extern UNKNOWN IMPORTED GLOBAL)
  if(WIN32)
    set_target_properties(feat-cudss-extern PROPERTIES IMPORTED_LOCATION ${cudss_SOURCE_DIR}/lib/cudss.lib)
  else()
    set_target_properties(feat-cudss-extern PROPERTIES IMPORTED_LOCATION ${cudss_SOURCE_DIR}/lib/libcudss.so)
  endif()
  target_include_directories(feat-cudss-extern INTERFACE ${cudss_SOURCE_DIR}/include)

  add_library(cuDSS::cuDSS ALIAS feat-cudss-extern)
endif()
