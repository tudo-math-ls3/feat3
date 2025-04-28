# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

set(kernel-cubature-list
  empty_cubature.cpp
)

if (FEAT_EICKT)
  set(kernel-cubature-list
    ${kernel-cubature-list}
    dynamic_factory-eickt.cpp
  )
endif (FEAT_EICKT)

add_library(feat-kernel-cubature STATIC ${kernel-cubature-list})

target_link_libraries(feat-kernel-cubature PRIVATE feat-kernel-util)

# list of cubature tests
set(test_list
)

include(FEATKernelTests)
add_feat_kernel_tests(TESTS test_list TARGET cubature_tests)
