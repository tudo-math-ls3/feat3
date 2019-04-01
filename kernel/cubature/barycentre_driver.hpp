// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
#define KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/driver_base.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAT
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      class BarycentreDriverBase :
        public DriverBase<Shape_>
      {
      public:
        /// this rule is not variadic
        static constexpr bool variadic = false;
        /// this rule has one point
        static constexpr int num_points = 1;

        ///Returns the name of the cubature rule.
        static String name()
        {
          return "barycentre";
        }

        /**
         * \brief Adds the driver's aliases.
         *
         * \param[in] functor
         * The functor whose \p alias function is to be called.
         */
        template<typename Functor_>
        static void alias(Functor_& functor)
        {
          functor.alias("midpoint");
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Barycentre Rule driver class template
     *
     * This driver implements the barycentre cubature rule.
     * \see http://en.wikipedia.org/wiki/Rectangle_method
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class BarycentreDriver DOXY({});

    // Simplex specialisation
    template<int dim_>
    class BarycentreDriver<Shape::Simplex<dim_> > :
      public Intern::BarycentreDriverBase<Shape::Simplex<dim_> >
    {
    public:
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
        rule.get_weight(0) = Weight_(1) / Weight_(MetaMath::Factorial<dim_>::value);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = Coord_(1) / Coord_(dim_ + 1);
        }
      }
    }; // class BarycentreDriver<Simplex<...>>

    // Hypercube specialisation
    template<int dim_>
    class BarycentreDriver<Shape::Hypercube<dim_> > :
      public Intern::BarycentreDriverBase<Shape::Hypercube<dim_> >
    {
    public:
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
        rule.get_weight(0) = Weight_(1 << dim_);

        // create coords of barycentre point
        for(int i(0); i < dim_; ++i)
        {
          rule.get_coord(0, i) = Coord_(0);
        }
      }
    }; // class BarycentreDriver<Hypercube<...>>
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_BARYCENTRE_DRIVER_HPP
