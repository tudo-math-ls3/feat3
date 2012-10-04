#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Maclaurin Rule driver class template
       *
       * This driver implements the Maclaurin cubature rules.
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
      class MaclaurinDriver :
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
          /// this rule has at least 1 points
          min_points = 1,
          /// this rule has at most 5 points
          max_points = 5
        };

        ///Returns the name of the cubature rule.
        static String name()
        {
          return "maclaurin";
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
          case 1:
            create_1(rule);
            break;
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
          }
        }

        /// Create 1-point Maclaurin rule (alias barycentre rule)
        static void create_1(RuleType& rule)
        {
          rule.get_coord(0) = Coord_(0);

          rule.get_weight(0) = WeightType(2);
        }

        /// Create 2-point Maclaurin rule
        static void create_2(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1) / Coord_(2);
          rule.get_coord(1) =  Coord_(1) / Coord_(2);

          rule.get_weight(0) = WeightType(1);
          rule.get_weight(1) = WeightType(1);
        }

        /// Create 3-point Maclaurin rule
        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(2) / Coord_(3);
          rule.get_coord(1) = Coord_(0);
          rule.get_coord(2) = Coord_(2) / Coord_(3);

          rule.get_weight(0) = WeightType(3) / WeightType(4);
          rule.get_weight(1) = WeightType(1) / WeightType(2);
          rule.get_weight(2) = WeightType(3) / WeightType(4);
        }

        /// Create 4-point Maclaurin rule
        static void create_4(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(3) / Coord_(4);
          rule.get_coord(1) = -Coord_(1) / Coord_(4);
          rule.get_coord(2) =  Coord_(1) / Coord_(4);
          rule.get_coord(3) =  Coord_(3) / Coord_(4);

          rule.get_weight(0) = Weight_(13) / Weight_(24);
          rule.get_weight(1) = Weight_(11) / Weight_(24);
          rule.get_weight(2) = Weight_(11) / Weight_(24);
          rule.get_weight(3) = Weight_(13) / Weight_(24);
        }

        /// Create 5-point Maclaurin rule
        static void create_5(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(4) / Coord_(5);
          rule.get_coord(1) = -Coord_(2) / Coord_(5);
          rule.get_coord(2) =  Coord_(0);
          rule.get_coord(3) =  Coord_(2) / Coord_(5);
          rule.get_coord(4) =  Coord_(4) / Coord_(5);

          rule.get_weight(0) = Weight_(275) / Weight_(576);
          rule.get_weight(1) = Weight_(100) / Weight_(576);
          rule.get_weight(2) = Weight_(402) / Weight_(576);
          rule.get_weight(3) = Weight_(100) / Weight_(576);
          rule.get_weight(4) = Weight_(275) / Weight_(576);
        }

      }; // class MaclaurinDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP