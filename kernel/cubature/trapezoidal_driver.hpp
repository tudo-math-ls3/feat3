#pragma once
#ifndef KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP
#define KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/rule.hpp>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      class TrapezoidalDriverBase
      {
      public:
        enum
        {
          variadic = 0,
        };

        static String name()
        {
          return "trapezoidal";
        }
      };
    } // namespace Intern
    /// \endcond

    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TrapezoidalDriver;

    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TrapezoidalDriver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::TrapezoidalDriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = dim_ + 1
      };

      static void create(RuleType& rule)
      {
        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = Weight_(1) / Weight_(Factorial<dim_ + 1>::value);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = Index(j+1) == i ? Coord_(1) : Coord_(0);
          }
        }
      }
    }; // class TrapezoidalDriver<Simplex<...>,...>

    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class TrapezoidalDriver<Shape::Hypercube<dim_>, Weight_, Coord_, Point_> :
      public Intern::TrapezoidalDriverBase
    {
    public:
      typedef Rule<Shape::Hypercube<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = (1 << dim_)
      };

      static void create(RuleType& rule)
      {
        for(Index i(0); i < Index(1 << dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = Weight_(1);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = Coord_(((i >> j) & 1) << 1) - Coord_(1);
          }
        }
      }
    }; // class TrapezoidalDriver<Hypercube<...>,...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP
