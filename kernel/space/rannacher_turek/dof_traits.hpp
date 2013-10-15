#pragma once
#ifndef KERNEL_SPACE_RANNACHER_TUREK_DOF_TRAITS_HPP
#define KERNEL_SPACE_RANNACHER_TUREK_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/space/rannacher_turek/variant.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace RannacherTurek
    {
      template<typename Shape_, typename Variant_>
      struct DofTag
      {
      };

      /**
       * \brief Lagrange-1 Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Tag_, int dim_>
      struct DofTraits
      {
        enum
        {
          /// no dofs for any shape dimension > 0
          count = 0
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };

      template<>
      struct DofTraits<DofTag<Shape::Hypercube<2>, Variant::StdNonPar>, 1>
      {
        enum
        {
          count = 1
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };

      template<>
      struct DofTraits<DofTag<Shape::Hypercube<3>, Variant::StdNonPar>, 2>
      {
        enum
        {
          count = 1
        };

        static Index derive_order(Index)
        {
          return Index(0);
        }
      };
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_DOF_TRAITS_HPP
