// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/cubature/scalar/rule.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Scalar cubature driver base class.
       *
       * \author Peter Zajac
       */
      class DriverBase
      {
      public:
        /// by default, tensorize the cubature forumula
        static constexpr bool tensorize = true;

        /**
         * \brief Applies an alias-functor.
         */
        template<typename Functor_>
        static void alias(Functor_&)
        {
          // do nothing
        }
      }; // class DriverBase
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT
