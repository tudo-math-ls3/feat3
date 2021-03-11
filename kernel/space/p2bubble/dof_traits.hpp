// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_P2BUBBLE_DOF_TRAITS_HPP
#define KERNEL_SPACE_P2BUBBLE_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace P2Bubble
    {
      /**
       * \brief P2-Bubble Dof-Traits class template.
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

      template<int shape_dim_>
      struct DofTraits<Shape::Simplex<shape_dim_>, 1>
      {
        /// 1 dof per simplex edge
        static constexpr int count = 1;
      };

      template<>
      struct DofTraits<Shape::Simplex<2>, 2>
      {
        /// 1 dof per 2D triangle
        static constexpr int count = 1;
      };

      template<>
      struct DofTraits<Shape::Simplex<3>, 3>
      {
        /// 1 dof per 3D tetrahedron
        static constexpr int count = 1;
      };
    } // namespace P2Bubble
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_P2BUBBLE_DOF_TRAITS_HPP
