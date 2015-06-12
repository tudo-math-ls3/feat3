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
        /// no dofs for any shape dimension > 0
        static constexpr int count = 0;
      };

      template<>
      struct DofTraits<DofTag<Shape::Hypercube<2>, Variant::StdNonPar>, 1>
      {
        static constexpr int count = 1;
      };

      template<>
      struct DofTraits<DofTag<Shape::Hypercube<3>, Variant::StdNonPar>, 2>
      {
        static constexpr int count = 1;
      };
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_DOF_TRAITS_HPP
