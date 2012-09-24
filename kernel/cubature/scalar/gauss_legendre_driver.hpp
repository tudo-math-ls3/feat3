#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/rule.hpp>

// includes, STL
#include <cmath>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      template<
        typename Weight_,
        typename Coord_>
      class GaussLegendreDriver
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<WeightType, CoordType> RuleType;

        enum
        {
          variadic = 1,
          tensorise = 1,
          min_points = 1,
          max_points = 5
        };

        static String name()
        {
          return "gauss-legendre";
        }

        static void create(Index num_points, RuleType& rule)
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

        static void create_1(RuleType& rule)
        {
          rule.get_coord(0) = CoordType(0);
          rule.get_weight(0) = WeightType(2);
        }

        static void create_2(RuleType& rule)
        {
          rule.get_coord(0) = -std::sqrt(CoordType(1) / CoordType(3));
          rule.get_coord(1) = +std::sqrt(CoordType(1) / CoordType(3));

          rule.get_weight(0) = WeightType(1);
          rule.get_weight(1) = WeightType(1);
        }

        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = -std::sqrt(CoordType(3) / CoordType(5));
          rule.get_coord(1) = CoordType(0);
          rule.get_coord(2) = +std::sqrt(CoordType(3) / CoordType(5));

          rule.get_weight(0) = WeightType(5) / WeightType(9);
          rule.get_weight(1) = WeightType(8) / WeightType(9);
          rule.get_weight(2) = WeightType(5) / WeightType(9);
        }

        static void create_4(RuleType& rule)
        {
          const CoordType dc = std::sqrt(CoordType(24) / CoordType(5));
          rule.get_coord(0) = -std::sqrt((CoordType(3) + dc) / CoordType(7));
          rule.get_coord(1) = -std::sqrt((CoordType(3) - dc) / CoordType(7));
          rule.get_coord(2) = +std::sqrt((CoordType(3) - dc) / CoordType(7));
          rule.get_coord(3) = +std::sqrt((CoordType(3) + dc) / CoordType(7));

          const WeightType dw = std::sqrt(WeightType(30));
          rule.get_weight(0) = (WeightType(18) - dw) / WeightType(36);
          rule.get_weight(1) = (WeightType(18) + dw) / WeightType(36);
          rule.get_weight(2) = (WeightType(18) + dw) / WeightType(36);
          rule.get_weight(3) = (WeightType(18) - dw) / WeightType(36);
        }

        static void create_5(RuleType& rule)
        {
          const CoordType dc = CoordType(2) * std::sqrt(CoordType(10) / CoordType(7));
          rule.get_coord(0) = -std::sqrt(CoordType(5) + dc) / CoordType(3);
          rule.get_coord(1) = -std::sqrt(CoordType(5) - dc) / CoordType(3);
          rule.get_coord(2) = CoordType(0);
          rule.get_coord(3) = +std::sqrt(CoordType(5) - dc) / CoordType(3);
          rule.get_coord(4) = +std::sqrt(CoordType(5) + dc) / CoordType(3);

          const WeightType dw = WeightType(13) * std::sqrt(WeightType(70));
          rule.get_weight(0) = (WeightType(322) - dw) / WeightType(900);
          rule.get_weight(1) = (WeightType(322) + dw) / WeightType(900);
          rule.get_weight(2) =  WeightType(128)       / WeightType(225);
          rule.get_weight(3) = (WeightType(322) + dw) / WeightType(900);
          rule.get_weight(4) = (WeightType(322) - dw) / WeightType(900);
        }
      }; // class GaussLegendreDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
