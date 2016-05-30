#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Midpoint Rule driver class template
       *
       * This driver implements the midpoint rule.
       *
       * \author Peter Zajac
       */
      class MidpointDriver :
        public DriverBase
      {
      public:
        /// this rule is not variadic
        static constexpr bool variadic = false;
        /// disable tensorisation
        static constexpr bool tensorise = false;
        /// this rule has 1 point
        static constexpr int num_points = 1;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "midpoint";
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
          rule.get_coord(0) = Coord_(0);
          rule.get_weight(0) = Weight_(2);
        }
      }; // class MidpointDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP
