#pragma once
#ifndef KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
#define KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP 1

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
      class BarycentreDriverBase
      {
      public:
        enum
        {
          variadic = 0,
          num_points = 1
        };

        static String name()
        {
          return "barycentre";
        }
      };
    } // namespace Intern
    /// \endcond

    template<
      typename Policy_,
      typename Shape_ = typename Policy_::ShapeType>
    class BarycentreDriver;

    template<
      typename Policy_,
      int dim_>
    class BarycentreDriver<Policy_, Shape::Simplex<dim_> > :
      public Intern::BarycentreDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;

      static void create(RuleType& rule)
      {
        rule.get_weight(0) = WeightType(1) / WeightType(Factorial<dim_>::value);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = CoordType(1) / CoordType(dim_ + 1);
        }
      }
    }; // class BarycentreDriver<...,Simplex<...>>

    template<
      typename Policy_,
      int dim_>
    class BarycentreDriver<Policy_, Shape::Hypercube<dim_> > :
      public Intern::BarycentreDriverBase<Policy_>
    {
    public:
      typedef Rule<Policy_> RuleType;
      typedef typename RuleType::WeightType WeightType;
      typedef typename RuleType::CoordType CoordType;

      static void create(RuleType& rule)
      {
        rule.get_weight(0) = WeightType(1 << dim_);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = CoordType(0);
        }
      }
    }; // class BarycentreDriver<...,Hypercube<...>>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
