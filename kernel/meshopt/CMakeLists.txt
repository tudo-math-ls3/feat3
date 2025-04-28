# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

# list of meshopt tests
set(kernel-meshopt-list
    empty_meshopt.cpp
    )

if (FEAT_EICKT)
  set(kernel-meshopt-list
    ${kernel-meshopt-list}
      dudv_functional-eickt.cpp
      hyperelasticity_functional-eickt.cpp
      rumpf_functional-eickt.cpp
      rumpf_trafo-eickt.cpp
    )
endif (FEAT_EICKT)

add_library(feat-kernel-meshopt STATIC ${kernel-meshopt-list})

target_link_libraries(feat-kernel-meshopt
PUBLIC
  feat-dependency
PRIVATE
  feat-kernel-cubature feat-kernel-geometry feat-kernel-lafem feat-kernel-util)

# list of meshopt tests
set(test_list
    hyperelasticity_functional-test
    )

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET meshopt_tests)
