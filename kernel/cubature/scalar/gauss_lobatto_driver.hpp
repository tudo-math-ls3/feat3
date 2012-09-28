#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_base.hpp>

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
      class GaussLobattoDriver :
        public DriverBase
      {
      public:
        typedef Weight_ WeightType;
        typedef Coord_ CoordType;
        typedef Rule<WeightType, CoordType> RuleType;

        enum
        {
          variadic = 1,
          min_points = 3,
          max_points = 6
        };

        static String name()
        {
          return "gauss-lobatto";
        }

        static void fill(RuleType& rule, Index num_points)
        {
          // create the formula
          switch(num_points)
          {
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
          }
        }

        static void create_3(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = CoordType(0);
          rule.get_coord(2) = CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(3);
          rule.get_weight(1) = WeightType(4) / WeightType(3);
          rule.get_weight(2) = WeightType(1) / WeightType(3);
        }

        static void create_4(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt(CoordType(5)) / CoordType(5);
          rule.get_coord(2) = +std::sqrt(CoordType(5)) / CoordType(5);
          rule.get_coord(3) =  CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(6);
          rule.get_weight(1) = WeightType(5) / WeightType(6);
          rule.get_weight(2) = WeightType(5) / WeightType(6);
          rule.get_weight(3) = WeightType(1) / WeightType(6);
        }

        static void create_5(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt(CoordType(21)) / CoordType(7);
          rule.get_coord(2) = CoordType(0);
          rule.get_coord(3) = +std::sqrt(CoordType(21)) / CoordType(7);
          rule.get_coord(4) = CoordType(1);

          rule.get_weight(0) = WeightType(1) / WeightType(10);
          rule.get_weight(1) = WeightType(49) / WeightType(90);
          rule.get_weight(2) = WeightType(32) / WeightType(45);
          rule.get_weight(3) = WeightType(49) / WeightType(90);
          rule.get_weight(4) = WeightType(1) / WeightType(10);
        }

        static void create_6(RuleType& rule)
        {
          rule.get_coord(0) = - CoordType(1);
          rule.get_coord(1) = -std::sqrt((CoordType(7) + std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(2) = -std::sqrt((CoordType(7) - std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(3) = std::sqrt((CoordType(7) - std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(4) = std::sqrt((CoordType(7) + std::sqrt(CoordType(28)))/ CoordType(21));
          rule.get_coord(5) = CoordType(1);


          rule.get_weight(0) = WeightType(1) / WeightType(15);
          rule.get_weight(1) = (WeightType(14) - std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(2) = (WeightType(14) + std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(3) = (WeightType(14) + std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(4) = (WeightType(14) - std::sqrt(WeightType(7))) / WeightType(30);
          rule.get_weight(5) = WeightType(1) / WeightType(15);
        }

      }; // class GaussLobattoDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP
