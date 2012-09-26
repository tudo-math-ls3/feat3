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
      /**
       * \brief Pulcherrima Rule driver class template
       *
       * This driver implements the Puilcherrima rule, a.k.a. 3/8-Simpson rule, see
       * \see http://en.wikipedia.org/wiki/Simpson's_rule
       * \see http://mathworld.wolfram.com/Simpsons38Rule.html
       *
       * This rule is also known as:
       * - 3rd order closed Newton-Cotes formula, see
       *   http://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas
       * - 3/8-Simpson rule
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
      class PulcherrimaDriver :
        public DriverBase
      {
      public:
        /// dummy enum
        enum
        {
          /// this rule is not variadic
          variadic = 0,
          /// this rule has 4 points
          num_points = 4
        };

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "pulcherrima";
        }

        /**
         * \brief Adds the driver's aliases.
         *
         * \param[in] functor
         * The functor whose \p alias function is to be called.
         */
        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          // insert all aliases for this cubature rule
          functor.alias("simpson-3/8");
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         */
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
