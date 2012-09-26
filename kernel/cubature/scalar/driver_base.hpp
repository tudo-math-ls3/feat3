#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_DRIVER_BASE_HPP
#define KERNEL_CUBATURE_SCALAR_DRIVER_BASE_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      class DriverBase
      {
      public:
        enum
        {
          tensorise = 1
        };

        template<typename Functor_>
        static void alias(Functor_&)
        {
          // do nothing
        }
      };
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_BASE_HPP