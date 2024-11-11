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
       * \brief Midpoint Rule driver class template
       *
       * This driver implements the midpoint rule.
       *
       * \author Peter Zajac
       */
      class MidpointDriver :
        public DriverBase
      {
      public:
        /// this rule is not variadic
        static constexpr bool variadic = false;
        /// disable tensorization
        static constexpr bool tensorize = false;
        /// this rule has 1 point
        static constexpr int num_points = 1;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "midpoint";
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
          rule.get_coord(0) = Coord_(0);
          rule.get_weight(0) = Weight_(2);
        }
      }; // class MidpointDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT
