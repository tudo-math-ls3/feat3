#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP 1

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
      class GaussLegendreDriver :
        public DriverBase
      {
      public:
        enum
        {
          variadic = 1,
          min_points = 1,
          max_points = 5
        };

        static String name()
        {
          return "gauss-legendre";
        }

        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          // midpoint rule -> gauss-legendre:1
          functor.alias("midpoint", 1);
        }

        static void fill(Rule<Weight_, Coord_>& rule, Index num_points)
        {
          // auxiliary variables
          Coord_ dc;
          Weight_ dw;

          // how many points do we have?
          switch(num_points)
          {
          case 1:
            rule.get_coord(0) = Coord_(0);

            rule.get_weight(0) = Weight_(2);
            break;

          case 2:
            rule.get_coord(0) = -std::sqrt(Coord_(1) / Coord_(3));
            rule.get_coord(1) = +std::sqrt(Coord_(1) / Coord_(3));

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -std::sqrt(Coord_(3) / Coord_(5));
            rule.get_coord(1) = Coord_(0);
            rule.get_coord(2) = +std::sqrt(Coord_(3) / Coord_(5));

            rule.get_weight(0) = Weight_(5) / Weight_(9);
            rule.get_weight(1) = Weight_(8) / Weight_(9);
            rule.get_weight(2) = Weight_(5) / Weight_(9);
            break;

          case 4:
            dc = std::sqrt(Coord_(24) / Coord_(5));
            rule.get_coord(0) = -std::sqrt((Coord_(3) + dc) / Coord_(7));
            rule.get_coord(1) = -std::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(2) = +std::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(3) = +std::sqrt((Coord_(3) + dc) / Coord_(7));

            dw = std::sqrt(Weight_(30));
            rule.get_weight(0) = (Weight_(18) - dw) / Weight_(36);
            rule.get_weight(1) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(2) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(3) = (Weight_(18) - dw) / Weight_(36);
            break;

          case 5:
            dc = Coord_(2) * std::sqrt(Coord_(10) / Coord_(7));
            rule.get_coord(0) = -std::sqrt(Coord_(5) + dc) / Coord_(3);
            rule.get_coord(1) = -std::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(2) = Coord_(0);
            rule.get_coord(3) = +std::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(4) = +std::sqrt(Coord_(5) + dc) / Coord_(3);

            dw = Weight_(13) * std::sqrt(Weight_(70));
            rule.get_weight(0) = (Weight_(322) - dw) / Weight_(900);
            rule.get_weight(1) = (Weight_(322) + dw) / Weight_(900);
            rule.get_weight(2) =  Weight_(128)       / Weight_(225);
            rule.get_weight(3) = (Weight_(322) + dw) / Weight_(900);
            rule.get_weight(4) = (Weight_(322) - dw) / Weight_(900);
            break;
          }
        }
      }; // class GaussLegendreDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
