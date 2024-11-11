// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Trapezoidal Rule driver class template
       *
       * This driver implements the trapezoidal rule.
       * \see http://en.wikipedia.org/wiki/Trapezoidal_rule
       * \see http://mathworld.wolfram.com/TrapezoidalRule.html
       *
       * \author Peter Zajac
       */
      class TrapezoidalDriver :
        public DriverBase
      {
      public:
        /// this rule is not variadic
        static constexpr bool variadic = false;
        /// disable tensorization
        static constexpr bool tensorize = false;
        /// this rule has 2 points
        static constexpr int num_points = 2;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "trapezoidal";
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         */
        template<
          typename Weight_,
          typename Coord_>
        static void fill(Rule<Weight_, Coord_>& rule)
        {
          rule.get_coord(0) = -Coord_(1);
          rule.get_coord(1) =  Coord_(1);
          rule.get_weight(0) = Weight_(1);
          rule.get_weight(1) = Weight_(1);
        }
      }; // class TrapezoidalDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT
