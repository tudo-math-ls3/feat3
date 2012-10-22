#pragma once
#ifndef KERNEL_FOUNDATION_BASE_HH
#define KERNEL_FOUNDATION_BASE_HH 1

#include<kernel/base_header.hpp>

namespace FEAST
{
  namespace Foundation
  {
    ///polytope level identifiers
    enum PolytopeLevels
    {
      pl_vertex = 0,
      pl_edge,
      pl_face,
      pl_polyhedron
    };

    enum Dimensions
    {
      dim_1D = 1,
      dim_2D = 2,
      dim_3D = 3
    };
  }
}
#endif
