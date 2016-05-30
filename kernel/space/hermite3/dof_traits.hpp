#pragma once
#ifndef KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP
#define KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Hermite3
    {
      /**
       * \brief Hermite-3 Dof-Traits class template.
       *
       * \author Peter Zajac
       */
      template<typename Shape_, int dim_>
      struct DofTraits
      {
        static constexpr int count = 0;
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Simplex<shape_dim_>, 0>
      {
        /// n+1 dofs per vertex
        static constexpr int count = shape_dim_ + 1;
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Simplex<shape_dim_>, shape_dim_>
      {
        /// (n+3 over 3) - (n+1)^2 dofs per cell
        static constexpr int count = MetaMath::Binomial<shape_dim_ + 3, 3>::value - (shape_dim_ + 1)*(shape_dim_ + 1);
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Hypercube<shape_dim_>, 0>
      {
        /// n+1 dofs per vertex
        static constexpr int count = shape_dim_ + 1;
      };

      template<int shape_dim_>
      struct DofTraits<Shape::Hypercube<shape_dim_>, shape_dim_>
      {
        /// (2^n - (n+1))*2^n dofs per cell
        static constexpr int count = ((1 << shape_dim_) - (shape_dim_ + 1)) * (1 << shape_dim_);
      };
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_HERMITE3_DOF_TRAITS_HPP
