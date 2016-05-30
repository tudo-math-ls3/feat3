#pragma once
#ifndef KERNEL_SPACE_LAGRANGE1_DOF_TRAITS_HPP
#define KERNEL_SPACE_LAGRANGE1_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Lagrange1
    {
      /**
       * \brief Lagrange-1 Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        /// no dofs for any shape dimension > 0
        static constexpr int count = 0;
      };

      template<typename Shape_>
      struct DofTraits<Shape_, 0>
      {
        /// 1 dof per vertex
        static constexpr int count = 1;
      };
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE1_DOF_TRAITS_HPP
