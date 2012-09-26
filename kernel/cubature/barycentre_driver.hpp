#pragma once
#ifndef KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
#define KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_base.hpp>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      class BarycentreDriverBase :
        public DriverBase
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
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class BarycentreDriver;

    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class BarycentreDriver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::BarycentreDriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;

      static void fill(RuleType& rule)
      {
        rule.get_weight(0) = Weight_(1) / Weight_(Factorial<dim_>::value);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = Coord_(1) / Coord_(dim_ + 1);
        }
      }
    }; // class BarycentreDriver<Simplex<...>,...>

    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class BarycentreDriver<Shape::Hypercube<dim_>, Weight_, Coord_, Point_> :
      public Intern::BarycentreDriverBase
    {
    public:
      typedef Rule<Shape::Hypercube<dim_>, Weight_, Coord_, Point_> RuleType;

      static void fill(RuleType& rule)
      {
        rule.get_weight(0) = Weight_(1 << dim_);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = Coord_(0);
        }
      }
    }; // class BarycentreDriver<Hypercube<...>,...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
