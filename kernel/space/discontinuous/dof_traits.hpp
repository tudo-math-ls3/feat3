#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_DOF_TRAITS_HPP
#define KERNEL_SPACE_DISCONTINUOUS_DOF_TRAITS_HPP 1

// includes, FEAST
#include <kernel/space/discontinuous/variant.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Discontinuous
    {
      template<typename Shape_, typename Variant_>
      struct DofTag
      {
      };

      /**
       * \brief Discontinuous Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Tag_, int dim_>
      struct DofTraits
      {
        static constexpr int count = 0;
      };

      template<int dim_, int degree_>
      struct DofTraits<DofTag<Shape::Hypercube<dim_>, Variant::StdPolyP<degree_> >, dim_>
      {
        static constexpr int count = MetaMath::Binomial<dim_ + degree_, degree_>::value;
      };

      template<int dim_, int degree_>
      struct DofTraits<DofTag<Shape::Simplex<dim_>, Variant::StdPolyP<degree_> >, dim_>
      {
        static constexpr int count = MetaMath::Binomial<dim_ + degree_, degree_>::value;
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_DOF_TRAITS_HPP
