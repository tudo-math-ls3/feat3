// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_DRIVER_BASE_HPP
#define KERNEL_CUBATURE_DRIVER_BASE_HPP 1

// includes, FEAT
#include <kernel/cubature/rule.hpp>

namespace FEAT
{
  namespace Cubature
  {
    template<typename Shape_>
    class DriverBase
    {
    public:
      template<typename Functor_>
      static void alias(Functor_&)
      {
        // do nothing
      }
    };

    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      class DummyDriver :
        public DriverBase<Shape_>
      {
      public:
        static constexpr bool variadic = false;
        static constexpr int num_points = 0;

        static String name()
        {
          return "dummy";
        }

        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          functor.alias("dummy-alias");
        }

        template<
          typename Weight_,
          typename Coord_,
          typename Point_>
        static void fill(Rule<Shape_, Weight_, Coord_, Point_>&)
        {
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_DRIVER_BASE_HPP
