// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP
#define KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP 1

// includes, FEAT
#include <kernel/space/base.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Discontinuous
    {
      /**
       * \brief Discontinuous element variant namespace
       */
      namespace Variant
      {
        /**
         * \brief Standard discontinous parametric P_k variant
         */
        template<int degree_ = 0>
        struct StdPolyP
        {
          /// This variant's local polynomial degree
          static constexpr int local_degree = degree_;

          static String name()
          {
            return String("StdPolyP<") + stringify(degree_) + ">";
          }
        };
      } // namespace Variant
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP
