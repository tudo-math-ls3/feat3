#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP
#define KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
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
          static String name()
          {
            return String("StdPolyP<") + stringify(degree_) + ">";
          }
        };
      } // namespace Variant
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_VARIANT_HPP
