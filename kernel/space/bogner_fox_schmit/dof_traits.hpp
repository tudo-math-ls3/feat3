#pragma once
#ifndef KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP
#define KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace BognerFoxSchmit
    {
      /**
       * \brief Bogner-Fox-Schmit Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        /// no dofs for any shape dimension > 0
        static constexpr int count = 0;
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Hypercube<shape_dim_>, 0>
      {
        /// 2^n dofs per vertex
        static constexpr int count = (1 << shape_dim_);
      };
    } // namespace BognerFoxSchmit
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP
