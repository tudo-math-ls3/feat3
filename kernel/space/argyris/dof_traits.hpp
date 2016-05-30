#pragma once
#ifndef KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP
#define KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
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
        static constexpr int count = 0;
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 0>
      {
        /// 6 dofs per vertex
        static constexpr int count = 6;
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 1>
      {
        /// 1 dof per edge
        static constexpr int count = 1;
      };
    } // namespace Argyris
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_ARGYRIS_DOF_TRAITS_HPP
