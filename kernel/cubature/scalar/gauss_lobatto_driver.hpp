// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP 1

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
       * \brief Gauss-Lobatto Rule driver class template
       *
       * This driver implements the Gauss-Lobatto rule.
       * \see http://en.wikipedia.org/wiki/Gaussian_quadrature
       *
       * \author Constantin Christof
       */
      class GaussLobattoDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 3 points
        static constexpr int min_points = 3;
        /// this rule has at most 6 points
        static constexpr int max_points = 6;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "gauss-lobatto";
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
          case 3:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) =  Coord_(0);
            rule.get_coord(2) =  Coord_(1);

            rule.get_weight(0) = Weight_(1) / Weight_(3);
            rule.get_weight(1) = Weight_(4) / Weight_(3);
            rule.get_weight(2) = Weight_(1) / Weight_(3);
            break;

          case 4:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Math::sqrt(Coord_(5)) / Coord_(5);
            rule.get_coord(2) = +Math::sqrt(Coord_(5)) / Coord_(5);
            rule.get_coord(3) = +Coord_(1);

            rule.get_weight(0) = Weight_(1) / Weight_(6);
            rule.get_weight(1) = Weight_(5) / Weight_(6);
            rule.get_weight(2) = Weight_(5) / Weight_(6);
            rule.get_weight(3) = Weight_(1) / Weight_(6);
            break;

          case 5:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Math::sqrt(Coord_(21)) / Coord_(7);
            rule.get_coord(2) =  Coord_(0);
            rule.get_coord(3) = +Math::sqrt(Coord_(21)) / Coord_(7);
            rule.get_coord(4) = +Coord_(1);

            rule.get_weight(0) = Weight_(1) / Weight_(10);
            rule.get_weight(1) = Weight_(49) / Weight_(90);
            rule.get_weight(2) = Weight_(32) / Weight_(45);
            rule.get_weight(3) = Weight_(49) / Weight_(90);
            rule.get_weight(4) = Weight_(1) / Weight_(10);
            break;

          case 6:
            rule.get_coord(0) = -Coord_(1);
            rule.get_coord(1) = -Math::sqrt((Coord_(7) + Math::sqrt(Coord_(28)))/ Coord_(21));
            rule.get_coord(2) = -Math::sqrt((Coord_(7) - Math::sqrt(Coord_(28)))/ Coord_(21));
            rule.get_coord(3) =  Math::sqrt((Coord_(7) - Math::sqrt(Coord_(28)))/ Coord_(21));
            rule.get_coord(4) =  Math::sqrt((Coord_(7) + Math::sqrt(Coord_(28)))/ Coord_(21));
            rule.get_coord(5) =  Coord_(1);


            rule.get_weight(0) = Weight_(1) / Weight_(15);
            rule.get_weight(1) = (Weight_(14) - Math::sqrt(Weight_(7))) / Weight_(30);
            rule.get_weight(2) = (Weight_(14) + Math::sqrt(Weight_(7))) / Weight_(30);
            rule.get_weight(3) = (Weight_(14) + Math::sqrt(Weight_(7))) / Weight_(30);
            rule.get_weight(4) = (Weight_(14) - Math::sqrt(Weight_(7))) / Weight_(30);
            rule.get_weight(5) = Weight_(1) / Weight_(15);
            break;
          }
        }
      }; // class GaussLobattoDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LOBATTO_DRIVER_HPP
