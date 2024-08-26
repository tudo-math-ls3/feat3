// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_Q1TBNP_DOF_TRAITS_HPP
#define KERNEL_SPACE_Q1TBNP_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Q1TBNP
    {
      /**
      * \brief Q1~bubble non-parametric Dof-Traits class template.
      *
      * \author Peter Zajac
      */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        /// no dofs for any shape dimension > 0
        static constexpr int count = 0;
      };

      template<int dim_>
      struct DofTraits<Shape::Hypercube<dim_+1>, dim_>
      {
        /// 1 dof per facet
        static constexpr int count = 1;
      };

      template<>
      struct DofTraits<Shape::Hypercube<2>, 2>
      {
        /// 1 bubble DOF in 2D
        static constexpr int count = 1;
      };

      template<>
      struct DofTraits<Shape::Hypercube<3>, 3>
      {
        /// 3 bubble DOFs in 3D
        static constexpr int count = 3;
      };
    } // namespace Q1TBNP
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_Q1TBNP_DOF_TRAITS_HPP
