#pragma once
#ifndef KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP
#define KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_base.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      class TrapezoidalDriverBase :
        public DriverBase<Shape_>
      {
      public:
        enum
        {
          /// this rule is not variadic
          variadic = 0,
        };

        ///Returns the name of the cubature rule.
        static String name()
        {
          return "trapezoidal";
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Trapezoidal driver class template
     *
     * This driver implements the trapezoidal rule.
     * \see http://en.wikipedia.org/wiki/Trapezoidal_rule
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class TrapezoidalDriver DOXY({});

    // Simplex specialisation
    template<int dim_>
    class TrapezoidalDriver<Shape::Simplex<dim_> > :
      public Intern::TrapezoidalDriverBase<Shape::Simplex<dim_> >
    {
    public:
      enum
      {
        num_points = dim_ + 1
      };

      /**
       * \brief Fills the cubature rule structure.
       *
       * \param[in,out] rule
       * The cubature rule to be filled.
       */
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_>& rule)
      {
        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = Weight_(1) / Weight_(MetaMath::Factorial<dim_ + 1>::value);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = Index(j+1) == i ? Coord_(1) : Coord_(0);
          }
        }
      }
    }; // class TrapezoidalDriver<Simplex<...>>

    // Hypercube specialisation
    template<int dim_>
    class TrapezoidalDriver<Shape::Hypercube<dim_> > :
      public Intern::TrapezoidalDriverBase<Shape::Hypercube<dim_> >
    {
    public:
      enum
      {
        num_points = (1 << dim_) // = 2^dim_
      };

      /**
       * \brief Fills the cubature rule structure.
       *
       * \param[in,out] rule
       * The cubature rule to be filled.
       */
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(Rule<Shape::Hypercube<dim_>, Weight_, Coord_, Point_>& rule)
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
    }; // class TrapezoidalDriver<Hypercube<...>>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_TRAPEZOIDAL_DRIVER_HPP
