// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP
#define KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
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
} // namespace FEAT

#endif // KERNEL_SPACE_BOGNER_FOX_SCHMIT_DOF_TRAITS_HPP
