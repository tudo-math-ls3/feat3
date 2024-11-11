// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/discontinuous/variant.hpp>

namespace FEAT
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
} // namespace FEAT
