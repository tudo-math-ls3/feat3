#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAST
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
       * \tparam Weight_
       * The data type for the cubature weights.
       *
       * \tparam Coord_
       * The data type for the cubature point coordinates.
       *
       * \author Constantin Christof
       */
      template<
        typename Weight_,
        typename Coord_>
      class NewtonCotesClosedDriver :
        public DriverBase
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<WeightType, CoordType> RuleType;

        enum
        {
          /// this rule is variadic
          variadic = 1,
          /// this rule has at least 2 points
          min_points = 2,
          /// this rule has at most 7 points
          max_points = 7
        };

        ///Returns the name of the cubature rule.
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
         * \param[in] num_points
         * The number of quadrature points.
         */
        static void fill(RuleType& rule, Index num_points)
        {
          // how many points do we have?
          switch(num_points)
          {
          case 2:
            create_2(rule);
            break;
          case 3:
            create_3(rule);
            break;
          case 4:
            create_4(rule);
            break;
          case 5:
            create_5(rule);
            break;
          case 6:
            create_6(rule);
            break;
          case 7:
            create_7(rule);
            break;
          }
        }

        /// Create 2-point closed Newton-Cotes rule (alias trapezoidal rule)
        static void create_2(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) = Coord_(1);

          rule.get_weight(0) = WeightType(1);
          rule.get_weight(1) = WeightType(1);
        }

        // Create 3-point closed Newton-Cotes rule (alias Simpson rule)
        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) = Coord_(0);
          rule.get_coord(2) = Coord_(1);

          rule.get_weight(0) = WeightType(1) / WeightType(3);
          rule.get_weight(1) = WeightType(2) / WeightType(3);
          rule.get_weight(2) = WeightType(1) / WeightType(3);
        }

        // Create 4-point closed Newton-Cotes rule (alias pulcherrima rule)
        static void create_4(RuleType& rule)
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

        /// Create 5-point closed Newton-Cotes rule (alias Milne-Bool rule)
        static void create_5(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) = -Coord_(1) / Coord_(2);
          rule.get_coord(2) =  Coord_(0);
          rule.get_coord(3) =  Coord_(1) / Coord_(2);
          rule.get_coord(4) =  Coord_(1);

          rule.get_weight(0) = Weight_(7) / Weight_(45);
          rule.get_weight(1) = Weight_(32) / Weight_(45);
          rule.get_weight(2) = Weight_(12) / Weight_(45);
          rule.get_weight(3) = Weight_(32) / Weight_(45);
          rule.get_weight(4) = Weight_(7) / Weight_(45);
        }

        /// Create 6-point closed Newton-Cotes rule
        static void create_6(RuleType& rule)
        {
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
        }

        /// Create 7-point closed Newton-Cotes rule (alias Weddle rule)
        static void create_7(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) = -Coord_(2) / Coord_(3);
          rule.get_coord(2) = -Coord_(1) / Coord_(3);
          rule.get_coord(3) =  Coord_(0);
          rule.get_coord(4) =  Coord_(1) / Coord_(3);
          rule.get_coord(5) =  Coord_(2) / Coord_(3);
          rule.get_coord(6) =  Coord_(1);

          rule.get_weight(0) = Weight_(41) / Weight_(420);
          rule.get_weight(1) = Weight_(216) / Weight_(420);
          rule.get_weight(2) = Weight_(27) / Weight_(420);
          rule.get_weight(3) = Weight_(272) / Weight_(420);
          rule.get_weight(4) = Weight_(27) / Weight_(420);
          rule.get_weight(5) = Weight_(216) / Weight_(420);
          rule.get_weight(6) = Weight_(41) / Weight_(420);
        }

      }; // class NewtonCotesClosedDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_NEWTON_COTES_CLOSED_DRIVER_HPP