if(NOT TARGET CGAL::CGAL)
  FetchContent_GetProperties(cgal)
  find_package(Boost 1.81 REQUIRED)

  add_library(feat-cgal-extern INTERFACE)
  target_include_directories(feat-cgal-extern INTERFACE "${cgal_SOURCE_DIR}/include")
  target_link_libraries(feat-cgal-extern INTERFACE Boost::boost)

  add_library(CGAL::CGAL ALIAS feat-cgal-extern)
endif()
