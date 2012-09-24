#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_TRAPEZOIDAL_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_TRAPEZOIDAL_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      template<
        typename Weight_,
        typename Coord_>
      class TrapezoidalDriver
      {
      public:
        enum
        {
          variadic = 0,
          tensorise = 0,
          num_points = 2
        };

        static String name()
        {
          return "trapezoidal";
        }

        static void create(Rule<Weight_, Coord_>& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) =  Coord_(1);
          rule.get_weight(0) = Weight_(1);
          rule.get_weight(1) = Weight_(1);
        }
      }; // class TrapezoidalDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_TRAPEZOIDAL_DRIVER_HPP
