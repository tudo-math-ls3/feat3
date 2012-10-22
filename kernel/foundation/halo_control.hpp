#pragma once
#ifndef KERNEL_FOUNDATION_HALO_CONTROL_HPP
#define KERNEL_FOUNDATION_HALO_CONTROL_HPP

#include <kernel/foundation/halo.hpp>
#include<kernel/geometry/conformal_mesh.hpp>

using namespace FEAST::Geometry;

namespace FEAST
{
  namespace Foundation
  {
    template<Dimensions dim_>
    struct HaloControl
    {
    };

    template<>
    struct HaloControl<dim_1D>
    {
    };

    template<>
    struct HaloControl<dim_2D>
    {
    };

    template<>
    struct HaloControl<dim_3D>
    {
    };
  }
}

#endif
