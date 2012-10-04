#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAST
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
       * \tparam Weight_
       * The data type for the cubature weights.
       *
       * \tparam Coord_
       * The data type for the cubature point coordinates.
       *
       * \author Peter Zajac
       */
      template<
        typename Weight_,
        typename Coord_>
      class MidpointDriver :
        public DriverBase
      {
      public:
        /// dummy enum
        enum
        {
          /// this rule is not variadic
          variadic = 0,
          /// this rule has 2 points
          num_points = 1,
          /// disable tensorisation
          tensorise = 0
        };

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
        static void fill(Rule<Weight_, Coord_>& rule)
        {
          rule.get_coord(0) = Coord_(0);
          rule.get_weight(0) = Weight_(2);
        }
      }; // class MidpointDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_MIDPOINT_DRIVER_HPP
