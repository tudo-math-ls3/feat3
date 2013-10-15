#pragma once
#ifndef KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP
#define KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Hermite3
    {
      /**
       * \brief Hermite-3 Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        enum
        {
          count = 0
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Simplex<shape_dim_>, 0>
      {
        enum
        {
          /// n+1 dofs per vertex
          count = shape_dim_ + 1
        };

        static Index derive_order(Index assign_idx)
        {
          return assign_idx == Index(0) ? Index(0) : Index(1);
        }
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Simplex<shape_dim_>, shape_dim_>
      {
        enum
        {
          /// (n+3 over 3) - (n+1)^2 dofs per cell
          count = MetaMath::Binomial<shape_dim_ + 3, 3>::value - (shape_dim_ + 1)*(shape_dim_ + 1)
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Hypercube<shape_dim_>, 0>
      {
        enum
        {
          /// n+1 dofs per vertex
          count = shape_dim_ + 1
        };

        static Index derive_order(Index assign_idx)
        {
          return assign_idx == Index(0) ? Index(0) : Index(1);
        }
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Hypercube<shape_dim_>, shape_dim_>
      {
        enum
        {
          /// (2^n - (n+1))*2^n dofs per cell
          count = ((1 << shape_dim_) - (shape_dim_ + 1)) * (1 << shape_dim_)
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP
