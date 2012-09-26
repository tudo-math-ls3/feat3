#pragma once
#ifndef KERNEL_CUBATURE_DRIVER_BASE_HPP
#define KERNEL_CUBATURE_DRIVER_BASE_HPP 1

// includes, FEAST
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
      class DriverBase
      {
      public:
        template<typename Functor_>
        static void alias(Functor_&)
        {
          // do nothing
        }
      };

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DRIVER_BASE_HPP
