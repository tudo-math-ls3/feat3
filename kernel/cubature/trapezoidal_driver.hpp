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
      template<typename Policy_>
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
      typename Policy_,
      typename Shape_ = typename Policy_::ShapeType>
    class TrapezoidalDriver;

    template<
      typename Policy_,
      int dim_>
    class TrapezoidalDriver<Policy_, Shape::Simplex<dim_> > :
      public Intern::TrapezoidalDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      enum
      {
        num_points = dim_ + 1
      };

      static void create(RuleType& rule)
      {
        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = WeightType(1) / WeightType(Factorial<dim_ + 1>::value);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = Index(j+1) == i ? CoordType(1) : CoordType(0);
          }
        }
      }
    }; // class TrapezoidalDriver<...,Simplex<...>>

    template<
      typename Policy_,
      int dim_>
    class TrapezoidalDriver<Policy_, Shape::Hypercube<dim_> > :
      public Intern::TrapezoidalDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;
      enum
      {
        num_points = (1 << dim_)
      };

      static void create(RuleType& rule)
      {
        for(Index i(0); i < Index(1 << dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = WeightType(1);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = CoordType(((i >> j) & 1) << 1) - CoordType(1);
          }
        }
      }
    }; // class TrapezoidalDriver<...,Hypercube<...>>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP
