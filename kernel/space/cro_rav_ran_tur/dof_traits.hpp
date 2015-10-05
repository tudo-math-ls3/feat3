#pragma once
#ifndef KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP
#define KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>

namespace FEAST
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
} // namespace FEAST

#endif // KERNEL_SPACE_CRO_RAV_RAN_TUR_DOF_TRAITS_HPP
