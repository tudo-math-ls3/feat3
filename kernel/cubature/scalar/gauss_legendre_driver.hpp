// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
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
       * \author Peter Zajac
       */
      class GaussLegendreDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 1 point
        static constexpr int min_points = 1;
        /// this rule has at most 5 points
        static constexpr int max_points = 5;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "gauss-legendre";
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
            rule.get_coord(0) = -Math::sqrt(Coord_(1) / Coord_(3));
            rule.get_coord(1) = +Math::sqrt(Coord_(1) / Coord_(3));

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -Math::sqrt(Coord_(3) / Coord_(5));
            rule.get_coord(1) = Coord_(0);
            rule.get_coord(2) = +Math::sqrt(Coord_(3) / Coord_(5));

            rule.get_weight(0) = Weight_(5) / Weight_(9);
            rule.get_weight(1) = Weight_(8) / Weight_(9);
            rule.get_weight(2) = Weight_(5) / Weight_(9);
            break;

          case 4:
            dc = Math::sqrt(Coord_(24) / Coord_(5));
            rule.get_coord(0) = -Math::sqrt((Coord_(3) + dc) / Coord_(7));
            rule.get_coord(1) = -Math::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(2) = +Math::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(3) = +Math::sqrt((Coord_(3) + dc) / Coord_(7));

            dw = Math::sqrt(Weight_(30));
            rule.get_weight(0) = (Weight_(18) - dw) / Weight_(36);
            rule.get_weight(1) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(2) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(3) = (Weight_(18) - dw) / Weight_(36);
            break;

          case 5:
            dc = Coord_(2) * Math::sqrt(Coord_(10) / Coord_(7));
            rule.get_coord(0) = -Math::sqrt(Coord_(5) + dc) / Coord_(3);
            rule.get_coord(1) = -Math::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(2) = Coord_(0);
            rule.get_coord(3) = +Math::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(4) = +Math::sqrt(Coord_(5) + dc) / Coord_(3);

            dw = Weight_(13) * Math::sqrt(Weight_(70));
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
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
