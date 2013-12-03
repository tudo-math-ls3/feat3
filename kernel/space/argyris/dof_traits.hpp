#pragma once
#ifndef KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP
#define KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Argyris
    {
      /**
       * \brief Argyris Dof-Traits class template.
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

      template<>
      struct DofTraits<Shape::Simplex<2>, 0>
      {
        enum
        {
          /// 6 dofs per vertex
          count = 6
        };

        static Index derive_order(Index assign_idx)
        {
          return assign_idx == Index(0) ? Index(0) : (assign_idx < 3  ? Index(1) : Index(2));
        }
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 1>
      {
        enum
        {
          /// 1 dof per edge
          count = 1
        };

        static Index derive_order(Index)
        {
          return Index(1);
        }
      };
    } // namespace Argyris
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP
