#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_INDEX_REFINER_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_INDEX_REFINER_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal/index_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Conformal
    {
      /// \cond internal
      namespace StandardRefinement
      {
        template<
          typename Shape_,
          int cell_dim_,
          int face_dim_>
        struct IndexRefiner;
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_INDEX_REFINER_HPP