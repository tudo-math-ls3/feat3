#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_TRAPEZOIDAL_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_TRAPEZOIDAL_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Trapezoidal Rule driver class template
       *
       * This driver implements the trapezoidal rule.
       * \see http://en.wikipedia.org/wiki/Trapezoidal_rule
       * \see http://mathworld.wolfram.com/TrapezoidalRule.html
       *
       * \author Peter Zajac
       */
      class TrapezoidalDriver :
        public DriverBase
      {
      public:
        /// this rule is not variadic
        static constexpr bool variadic = false;
        /// disable tensorisation
        static constexpr bool tensorise = false;
        /// this rule has 2 points
        static constexpr int num_points = 2;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "trapezoidal";
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         */
        template<
          typename Weight_,
          typename Coord_>
        static void fill(Rule<Weight_, Coord_>& rule)
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
