#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_PULCHERRIMA_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_PULCHERRIMA_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      template<
        typename Weight_,
        typename Coord_>
      class PulcherrimaDriver :
        public DriverBase
      {
      public:
        enum
        {
          variadic = 0,
          num_points = 4
        };

        static String name()
        {
          return "pulcherrima";
        }

        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          // insert all aliases for this cubature rule
          functor.alias("simpson-3/8");
        }

        static void fill(Rule<Weight_, Coord_>& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) = -Coord_(1) / Coord_(3);
          rule.get_coord(2) =  Coord_(1) / Coord_(3);
          rule.get_coord(3) =  Coord_(1);
          rule.get_weight(0) = Weight_(1) / Weight_(4);
          rule.get_weight(1) = Weight_(3) / Weight_(4);
          rule.get_weight(2) = Weight_(3) / Weight_(4);
          rule.get_weight(3) = Weight_(1) / Weight_(4);
        }
      }; // class PulcherrimaDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_PULCHERRIMA_DRIVER_HPP
