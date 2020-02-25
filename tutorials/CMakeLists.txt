# CmakeLists for directory "tutorials"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 3.9)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(tutorial-01-poisson tutorial_01_poisson.cpp)
TARGET_LINK_LIBRARIES(tutorial-01-poisson feat)

ADD_EXECUTABLE(tutorial-02-laplace tutorial_02_laplace.cpp)
TARGET_LINK_LIBRARIES(tutorial-02-laplace feat)

ADD_EXECUTABLE(tutorial-03-andicore tutorial_03_andicore.cpp)
TARGET_LINK_LIBRARIES(tutorial-03-andicore feat)

ADD_EXECUTABLE(tutorial-04-parser tutorial_04_parser.cpp)
TARGET_LINK_LIBRARIES(tutorial-04-parser feat)

ADD_EXECUTABLE(tutorial-05-multigrid tutorial_05_multigrid.cpp)
TARGET_LINK_LIBRARIES(tutorial-05-multigrid feat)

ADD_EXECUTABLE(tutorial-06-global tutorial_06_global.cpp)
TARGET_LINK_LIBRARIES(tutorial-06-global feat)

ADD_EXECUTABLE(tutorial-07-stokes tutorial_07_stokes.cpp)
TARGET_LINK_LIBRARIES(tutorial-07-stokes feat)

ADD_EXECUTABLE(tutorial-08-meshfile tutorial_08_meshfile.cpp)
TARGET_LINK_LIBRARIES(tutorial-08-meshfile feat)

# tutorials target
ADD_CUSTOM_TARGET(tutorials DEPENDS
  tutorial-01-poisson
  tutorial-02-laplace
  tutorial-03-andicore
  tutorial-04-parser
  tutorial-05-multigrid
  tutorial-06-global
  tutorial-07-stokes
  tutorial-08-meshfile
)
