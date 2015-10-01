# CmakeLists for directory "tutorials"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# applications to build
ADD_EXECUTABLE(tutorial-01-poisson tutorial_01_poisson.cpp)
TARGET_LINK_LIBRARIES(tutorial-01-poisson feast)

ADD_EXECUTABLE(tutorial-02-laplace tutorial_02_laplace.cpp)
TARGET_LINK_LIBRARIES(tutorial-02-laplace feast)

ADD_EXECUTABLE(tutorial-03-andicore tutorial_03_andicore.cpp)
TARGET_LINK_LIBRARIES(tutorial-03-andicore feast)

ADD_EXECUTABLE(tutorial-04-parser tutorial_04_parser.cpp)
TARGET_LINK_LIBRARIES(tutorial-04-parser feast)

ADD_EXECUTABLE(tutorial-x1-andicore-mg tutorial_x1_andicore_mg.cpp)
TARGET_LINK_LIBRARIES(tutorial-x1-andicore-mg feast)

# tutorials target
ADD_CUSTOM_TARGET(tutorials DEPENDS
  tutorial-01-poisson
  tutorial-02-laplace
  tutorial-03-andicore
  tutorial-04-parser
  tutorial-x1-andicore-mg
)
