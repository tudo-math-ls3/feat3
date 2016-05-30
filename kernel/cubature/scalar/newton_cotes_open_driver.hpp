#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_NEWTON_COTES_OPEN_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_NEWTON_COTES_OPEN_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Open Newton-Cotes Rule driver class template
       *
       * This driver implements the open Newton Cotes cubature rules.
       * \see http://de.wikipedia.org/wiki/Newton-Cotes-Formeln
       *
       * \author Constantin Christof
       */
      class NewtonCotesOpenDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 1 point
        static constexpr int min_points = 1;
        /// this rule has at most 7 points
        static constexpr int max_points = 7;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "newton-cotes-open";
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
          case 1:
            rule.get_coord(0) = Coord_(0);

            rule.get_weight(0) = Weight_(2);
            break;

          case 2:
            rule.get_coord(0) = -Coord_(1) / Coord_(3);
            rule.get_coord(1) =  Coord_(1) / Coord_(3);

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -Coord_(1) / Coord_(2);
            rule.get_coord(1) = Coord_(0);
            rule.get_coord(2) = Coord_(1) / Coord_(2);

            rule.get_weight(0) = Weight_(4) / Weight_(3);
            rule.get_weight(1) = - Weight_(2) / Weight_(3);
            rule.get_weight(2) = Weight_(4) / Weight_(3);
            break;

          case 4:
            rule.get_coord(0) = -Coord_(3) / Coord_(5);
            rule.get_coord(1) = -Coord_(1) / Coord_(5);
            rule.get_coord(2) =  Coord_(1) / Coord_(5);
            rule.get_coord(3) =  Coord_(3) / Coord_(5);

            rule.get_weight(0) = Weight_(11) / Weight_(12);
            rule.get_weight(1) = Weight_(1) / Weight_(12);
            rule.get_weight(2) = Weight_(1) / Weight_(12);
            rule.get_weight(3) = Weight_(11) / Weight_(12);
            break;

          case 5:
            rule.get_coord(0) = -Coord_(2) / Coord_(3);
            rule.get_coord(1) = -Coord_(1) / Coord_(3);
            rule.get_coord(2) =  Coord_(0);
            rule.get_coord(3) =  Coord_(1) / Coord_(3);
            rule.get_coord(4) =  Coord_(2) / Coord_(3);

            rule.get_weight(0) = Weight_(11) / Weight_(10);
            rule.get_weight(1) = - Weight_(14) / Weight_(10);
            rule.get_weight(2) = Weight_(26) / Weight_(10);
            rule.get_weight(3) = - Weight_(14) / Weight_(10);
            rule.get_weight(4) = Weight_(11) / Weight_(10);
            break;

          case 6:
            rule.get_coord(0) = -Coord_(5) / Coord_(7);
            rule.get_coord(1) = -Coord_(3) / Coord_(7);
            rule.get_coord(2) = -Coord_(1) / Coord_(7);
            rule.get_coord(3) =  Coord_(1) / Coord_(7);
            rule.get_coord(4) =  Coord_(3) / Coord_(7);
            rule.get_coord(5) =  Coord_(5) / Coord_(7);

            rule.get_weight(0) =  Weight_(611) / Weight_(720);
            rule.get_weight(1) = -Weight_(453) / Weight_(720);
            rule.get_weight(2) =  Weight_(562) / Weight_(720);
            rule.get_weight(3) =  Weight_(562) / Weight_(720);
            rule.get_weight(4) = -Weight_(453) / Weight_(720);
            rule.get_weight(5) =  Weight_(611) / Weight_(720);
            break;

          case 7:
            rule.get_coord(0) = -Coord_(3) / Coord_(4);
            rule.get_coord(1) = -Coord_(1) / Coord_(2);
            rule.get_coord(2) = -Coord_(1) / Coord_(4);
            rule.get_coord(3) =  Coord_(0);
            rule.get_coord(4) =  Coord_(1) / Coord_(4);
            rule.get_coord(5) =  Coord_(1) / Coord_(2);
            rule.get_coord(6) =  Coord_(3) / Coord_(4);

            rule.get_weight(0) =  Weight_(920)  / Weight_(945);
            rule.get_weight(1) = -Weight_(1908) / Weight_(945);
            rule.get_weight(2) =  Weight_(4392) / Weight_(945);
            rule.get_weight(3) = -Weight_(4918) / Weight_(945);
            rule.get_weight(4) =  Weight_(4392) / Weight_(945);
            rule.get_weight(5) = -Weight_(1908) / Weight_(945);
            rule.get_weight(6) =  Weight_(920)  / Weight_(945);
            break;
          }
        }
      }; // class NewtonCotesOpenDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_NEWTON_COTES_OPEN_DRIVER_HPP
