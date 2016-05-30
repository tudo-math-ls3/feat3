#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Closed Newton-Cotes Rule driver class template
       *
       * This driver implements the closed Newton Cotes cubature rules.
       * \see http://de.wikipedia.org/wiki/Newton-Cotes-Formeln
       *
       * \author Constantin Christof
       */
      class NewtonCotesClosedDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 2 points
        static constexpr int min_points = 2;
        /// this rule has at most 7 points
        static constexpr int max_points = 7;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "newton-cotes-closed";
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
          functor.alias("simpson", 3);
          functor.alias("pulcherrima", 4);
          functor.alias("milne-bool", 5);
          functor.alias("6-point", 6);
          functor.alias("weddle", 7);
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         *
         * \param[in] num_points
         * The number of quadrature points.
         */
        template<
          typename Weight_,
          typename Coord_>
        static void fill(Rule<Weight_, Coord_>& rule, int num_points)
        {
          // how many points do we have?
          switch(num_points)
          {
          case 2:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = Coord_(1);

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = Coord_(0);
            rule.get_coord(2) = Coord_(1);

            rule.get_weight(0) = Weight_(1) / Weight_(3);
            rule.get_weight(1) = Weight_(2) / Weight_(3);
            rule.get_weight(2) = Weight_(1) / Weight_(3);
            break;

          case 4:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Coord_(1) / Coord_(3);
            rule.get_coord(2) =  Coord_(1) / Coord_(3);
            rule.get_coord(3) =  Coord_(1);

            rule.get_weight(0) = Weight_(1) / Weight_(4);
            rule.get_weight(1) = Weight_(3) / Weight_(4);
            rule.get_weight(2) = Weight_(3) / Weight_(4);
            rule.get_weight(3) = Weight_(1) / Weight_(4);
            break;

          case 5:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Coord_(1) / Coord_(2);
            rule.get_coord(2) =  Coord_(0);
            rule.get_coord(3) =  Coord_(1) / Coord_(2);
            rule.get_coord(4) =  Coord_(1);

            rule.get_weight(0) = Weight_(7)  / Weight_(45);
            rule.get_weight(1) = Weight_(32) / Weight_(45);
            rule.get_weight(2) = Weight_(12) / Weight_(45);
            rule.get_weight(3) = Weight_(32) / Weight_(45);
            rule.get_weight(4) = Weight_(7)  / Weight_(45);
            break;

          case 6:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Coord_(3) / Coord_(5);
            rule.get_coord(2) = -Coord_(1) / Coord_(5);
            rule.get_coord(3) =  Coord_(1) / Coord_(5);
            rule.get_coord(4) =  Coord_(3) / Coord_(5);
            rule.get_coord(5) =  Coord_(1);

            rule.get_weight(0) = Weight_(19) / Weight_(144);
            rule.get_weight(1) = Weight_(75) / Weight_(144);
            rule.get_weight(2) = Weight_(50) / Weight_(144);
            rule.get_weight(3) = Weight_(50) / Weight_(144);
            rule.get_weight(4) = Weight_(75) / Weight_(144);
            rule.get_weight(5) = Weight_(19) / Weight_(144);
            break;

          case 7:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Coord_(2) / Coord_(3);
            rule.get_coord(2) = -Coord_(1) / Coord_(3);
            rule.get_coord(3) =  Coord_(0);
            rule.get_coord(4) =  Coord_(1) / Coord_(3);
            rule.get_coord(5) =  Coord_(2) / Coord_(3);
            rule.get_coord(6) =  Coord_(1);

            rule.get_weight(0) = Weight_(41)  / Weight_(420);
            rule.get_weight(1) = Weight_(216) / Weight_(420);
            rule.get_weight(2) = Weight_(27)  / Weight_(420);
            rule.get_weight(3) = Weight_(272) / Weight_(420);
            rule.get_weight(4) = Weight_(27)  / Weight_(420);
            rule.get_weight(5) = Weight_(216) / Weight_(420);
            rule.get_weight(6) = Weight_(41)  / Weight_(420);
            break;
          }
        }
      }; // class NewtonCotesClosedDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP
