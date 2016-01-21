# CmakeLists for directory "tutorials"

# this needs to be set again for each directory
cmake_minimum_required (VERSION 2.8)
#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

if (FEAST_HAVE_PARMETIS)
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/parmetis/parmetis-4.0.3/include/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/parmetis/parmetis-4.0.3/metis/include/")
endif (FEAST_HAVE_PARMETIS)

# applications to build
ADD_EXECUTABLE(tutorial-01-poisson tutorial_01_poisson.cpp)
TARGET_LINK_LIBRARIES(tutorial-01-poisson feast)

ADD_EXECUTABLE(tutorial-02-laplace tutorial_02_laplace.cpp)
TARGET_LINK_LIBRARIES(tutorial-02-laplace feast)

ADD_EXECUTABLE(tutorial-03-andicore tutorial_03_andicore.cpp)
TARGET_LINK_LIBRARIES(tutorial-03-andicore feast)

ADD_EXECUTABLE(tutorial-04-parser tutorial_04_parser.cpp)
TARGET_LINK_LIBRARIES(tutorial-04-parser feast)

ADD_EXECUTABLE(tutorial-05-multigrid tutorial_05_multigrid.cpp)
TARGET_LINK_LIBRARIES(tutorial-05-multigrid feast)

# tutorials target
ADD_CUSTOM_TARGET(tutorials DEPENDS
  tutorial-01-poisson
  tutorial-02-laplace
  tutorial-03-andicore
  tutorial-04-parser
  tutorial-05-multigrid
)
