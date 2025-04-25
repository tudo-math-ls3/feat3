if(NOT TARGET SuperLU::SuperLU)
  add_library(SuperLU::SuperLU ALIAS superlu_dist)
endif()
