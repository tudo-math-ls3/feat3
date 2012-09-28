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
          variadic = 1,
          min_points = 1,
          max_points = 5
        };

        static String name()
        {
          return "maclaurin";
        }

        static void fill(RuleType& rule, Index num_points)
        {
          // create the formula
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

        //midpoint
        static void create_1(RuleType& rule)
        {
          rule.get_coord(0) = Coord_(0);

          rule.get_weight(0) = WeightType(2);
        }

        static void create_2(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(1) / Coord_(2);
          rule.get_coord(1) =  Coord_(1) / Coord_(2);

          rule.get_weight(0) = WeightType(1);
          rule.get_weight(1) = WeightType(1);
        }

        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = -Coord_(2) / Coord_(3);
          rule.get_coord(1) = Coord_(0);
          rule.get_coord(2) = Coord_(2) / Coord_(3);

          rule.get_weight(0) = WeightType(3) / WeightType(4);
          rule.get_weight(1) = WeightType(1) / WeightType(2);
          rule.get_weight(2) = WeightType(3) / WeightType(4);
        }

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
