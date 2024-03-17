# CmakeLists for directory "tutorials"

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
add_executable(tutorial-01-poisson tutorial_01_poisson.cpp)
target_link_libraries(tutorial-01-poisson feat)

add_executable(tutorial-02-laplace tutorial_02_laplace.cpp)
target_link_libraries(tutorial-02-laplace feat)

add_executable(tutorial-03-andicore tutorial_03_andicore.cpp)
target_link_libraries(tutorial-03-andicore feat)

add_executable(tutorial-04-parser tutorial_04_parser.cpp)
target_link_libraries(tutorial-04-parser feat)

add_executable(tutorial-05-multigrid tutorial_05_multigrid.cpp)
target_link_libraries(tutorial-05-multigrid feat)

add_executable(tutorial-06-global tutorial_06_global.cpp)
target_link_libraries(tutorial-06-global feat)

add_executable(tutorial-07-stokes tutorial_07_stokes.cpp)
target_link_libraries(tutorial-07-stokes feat)

add_executable(tutorial-08-meshfile tutorial_08_meshfile.cpp)
target_link_libraries(tutorial-08-meshfile feat)

# tutorials target
add_custom_target(tutorials DEPENDS
  tutorial-01-poisson
  tutorial-02-laplace
  tutorial-03-andicore
  tutorial-04-parser
  tutorial-05-multigrid
  tutorial-06-global
  tutorial-07-stokes
  tutorial-08-meshfile
)
