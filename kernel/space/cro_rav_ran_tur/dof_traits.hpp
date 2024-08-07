// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP
#define KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace CroRavRanTur
    {
      /**
      * \brief Crouzeix-Raviart / Rannacher-Turek Dof-Traits class template.
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
      struct DofTraits<Shape::Simplex<dim_+1>, dim_>
      {
        /// 1 dof per facet
        static constexpr int count = 1;
      };

      template<int dim_>
      struct DofTraits<Shape::Hypercube<dim_+1>, dim_>
      {
        /// 1 dof per facet
        static constexpr int count = 1;
      };
    } // namespace CroRavRanTur
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP
