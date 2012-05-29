#pragma once
#ifndef KERNEL_GEOMETRY_CONF_STD_REF_TARGET_REFINER_HPP
#define KERNEL_GEOMETRY_CONF_STD_REF_TARGET_REFINER_HPP 1

// includes, FEAST
#include <kernel/geometry/target_set.hpp>
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
          int cell_dim_>
        struct TargetRefiner;

        template<>
        struct TargetRefiner<Shape::Vertex, 0>
        {
          typedef Shape::Vertex ShapeType;
          typedef TargetSet TargetSetType;
          typedef TargetSetHolder<ShapeType> TargetSetHolderType;

          static Index refine(
            TargetSetType& target_set_out,
            const Index offset,
            const Index* index_offsets,
            const TargetSetHolderType& target_set_holder_in)
          {
            const Index num_verts = target_set_holder_in.get_num_entities(0);
            const TargetSetType& target_set_in = target_set_holder_in.get_target_set<0>();

            for(Index i(0); i < num_verts; ++i)
            {
              target_set_out[i] = index_offsets[0] + target_set_in[offset + i];
            }

            return num_verts;;
          }
        };
      } // namespace StandardRefinement
      /// \endcond
    } // namespace Conformal
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONF_STD_REF_TARGET_REFINER_HPP
