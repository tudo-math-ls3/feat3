// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_DOF_TRAITS_HPP
#define KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace CaiDouSanSheYe
    {
      /**
      * \brief Cai-Douglas-Santos-Sheen-Ye Dof-Traits class template.
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

      template<int dim_>
      struct DofTraits<Shape::Hypercube<dim_>, dim_>
      {
        /// 1 dof per element
        static constexpr int count = 1;
      };
    } // namespace CaiDouSanSheYe
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_DOF_TRAITS_HPP
