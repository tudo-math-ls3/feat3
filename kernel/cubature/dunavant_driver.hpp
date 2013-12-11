#pragma once
#ifndef KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP
#define KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/symmetric_simplex_driver.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAST
{
  namespace Cubature
  {
    /**
     * \brief Silvester "open" driver class template
     *
     * This driver implements the open Dunavant rules for triangles.
     * \see D.A. Dunavant: High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle
     *      Int. j. numer. methods eng., Volume 21 (1985), pp. 1129 - 1148
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \todo implement rules for p = 3..13 and p > 14
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class DunavantDriver DOXY({});

    // Simplex<2> specialisation
    template<>
    class DunavantDriver<Shape::Simplex<2> > :
      public SymmetricSimplexDriver<Shape::Simplex<2> >
    {
    public:
      enum
      {
        /// this rule is variadic
        variadic = 1,
        min_points = 14,
        max_points = 14
      };

      /// Returns the name of the cubature rule.
      static String name()
      {
        return "dunavant";
      }

      static Index count(Index points)
      {
        switch(points)
        {
        case 14:
          return Index(42);
        case 15:
          return Index(48);
        default:
          return Index(0);
        }
      }

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
      static void fill(Rule<Shape::Simplex<2>, Weight_, Coord_, Point_>& rule, Index num_points)
      {
        Index off(0);
        switch(num_points)
        {
        case 14:
          off += fill_sym2(rule, off, Weight_(0.0109417906847145), Coord_(0.022072179275643), Coord_(0.488963910362179));
          off += fill_sym2(rule, off, Weight_(0.0163941767720625), Coord_(0.164710561319092), Coord_(0.417644719340454));
          off += fill_sym2(rule, off, Weight_(0.0258870522536460), Coord_(0.453044943382323), Coord_(0.273477528308839));
          off += fill_sym2(rule, off, Weight_(0.0210812943684965), Coord_(0.645588935174913), Coord_(0.177205532412543));
          off += fill_sym2(rule, off, Weight_(0.0072168498348885), Coord_(0.876400233818255), Coord_(0.061799883090873));
          off += fill_sym2(rule, off, Weight_(0.0024617018012000), Coord_(0.961218077502598), Coord_(0.019390961248701));
          off += fill_sym3(rule, off, Weight_(0.0123328766062820), Coord_(0.057124757403648), Coord_(0.172266687821356), Coord_(0.770608554774996));
          off += fill_sym3(rule, off, Weight_(0.0192857553935305), Coord_(0.092916249356972), Coord_(0.336861459796345), Coord_(0.570222290846683));
          off += fill_sym3(rule, off, Weight_(0.0072181540567670), Coord_(0.014646950055654), Coord_(0.298372882136258), Coord_(0.686980167808088));
          off += fill_sym3(rule, off, Weight_(0.0025051144192505), Coord_(0.001268330932872), Coord_(0.118974497696957), Coord_(0.879757171370171));
          break;
        }
      }
    }; // class DunavantDriver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP
