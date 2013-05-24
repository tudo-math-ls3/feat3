#pragma once
#ifndef KERNEL_CUBATURE_DRIVER_BASE_HPP
#define KERNEL_CUBATURE_DRIVER_BASE_HPP 1

// includes, FEAST
#include <kernel/cubature/rule.hpp>

namespace FEAST
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
        enum
        {
          variadic = 0,
          num_points = 0
        };

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
} // namespace FEAST

#endif // KERNEL_CUBATURE_DRIVER_BASE_HPP
