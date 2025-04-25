if(NOT TARGET deathhandler::deathhandler)
  FetchContent_GetProperties(death_handler)

  add_library(feat-deathhandler-extern STATIC ${death_handler_SOURCE_DIR}/death_handler.cc)
  target_include_directories(feat-deathhandler-extern INTERFACE SYSTEM ${death_handler_SOURCE_DIR})

  add_library(deathhandler::deathhandler ALIAS feat-deathhandler-extern)
endif()
