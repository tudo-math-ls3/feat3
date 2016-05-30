#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAT
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
       * \author Constantin Christof
       */
      class MaclaurinDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 1 points
        static constexpr int min_points = 1;
        /// this rule has at most 5 points
        static constexpr int max_points = 5;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "maclaurin";
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
            rule.get_coord(0) = -Coord_(1) / Coord_(2);
            rule.get_coord(1) =  Coord_(1) / Coord_(2);

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -Coord_(2) / Coord_(3);
            rule.get_coord(1) =  Coord_(0);
            rule.get_coord(2) =  Coord_(2) / Coord_(3);

            rule.get_weight(0) = Weight_(3) / Weight_(4);
            rule.get_weight(1) = Weight_(1) / Weight_(2);
            rule.get_weight(2) = Weight_(3) / Weight_(4);
            break;

          case 4:
            rule.get_coord(0) = -Coord_(3) / Coord_(4);
            rule.get_coord(1) = -Coord_(1) / Coord_(4);
            rule.get_coord(2) =  Coord_(1) / Coord_(4);
            rule.get_coord(3) =  Coord_(3) / Coord_(4);

            rule.get_weight(0) = Weight_(13) / Weight_(24);
            rule.get_weight(1) = Weight_(11) / Weight_(24);
            rule.get_weight(2) = Weight_(11) / Weight_(24);
            rule.get_weight(3) = Weight_(13) / Weight_(24);
            break;

          case 5:
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
            break;
          }
        }
      }; // class MaclaurinDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_MACLAURIN_DRIVER_HPP
