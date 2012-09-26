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
      /**
       * \brief Gauss-Legendre Rule driver class template
       *
       * This driver implements the Gauss-Legendre rule.
       * \see http://en.wikipedia.org/wiki/Gaussian_quadrature
       * \see http://mathworld.wolfram.com/GaussianQuadrature.html
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
      class GaussLegendreDriver :
        public DriverBase
      {
      public:
        /// dummy enum
        enum
        {
          /// this rule is variadic
          variadic = 1,
          /// this rule has at least 1 point
          min_points = 1,
          /// this rule has at most 5 points
          max_points = 5
        };

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "gauss-legendre";
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
          // midpoint rule -> gauss-legendre:1
          functor.alias("midpoint", 1);
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         */
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
