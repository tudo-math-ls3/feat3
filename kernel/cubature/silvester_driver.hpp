#pragma once
#ifndef KERNEL_CUBATURE_SILVESTER_DRIVER_HPP
#define KERNEL_CUBATURE_SILVESTER_DRIVER_HPP 1

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
     * This driver implements the open Silvester rules for simplices.
     * \see P. Silvester: Symmetric Quadrature Formulae for Simplexes;
     *      Math. Comp., Issue 24 (1970), pp. 95 - 100
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \todo implement 5,6 and 7 point rules for Simplex<2>
     * \todo implement Simplex<3> rules
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class SilvesterOpenDriver DOXY({});

    // Simplex<2> specialisation
    template<>
    class SilvesterOpenDriver<Shape::Simplex<2> > :
      public SymmetricSimplexDriver<Shape::Simplex<2> >
    {
    public:
      enum
      {
        /// this rule is variadic
        variadic = 1,
        min_points = 2,
        max_points = 8
      };

      /// Returns the name of the cubature rule.
      static String name()
      {
        return "silvester-open";
      }

      static Index count(Index points)
      {
        return ((points+1)*(points+2)) / Index(2);
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
        case 2:
          // 200
          off += fill_sym2(rule, off, Weight_( 7) / Weight_(24), Coord_(3) / Coord_(5), Coord_(1) / Coord_(5));
          // 011
          off += fill_sym2(rule, off, Weight_(-3) / Weight_(24), Coord_(1) / Coord_(5), Coord_(2) / Coord_(5));
          break;

        case 3:
          // 300
          off += fill_sym1(rule, off, Weight_(-12) / Weight_(60), Coord_(2) / Coord_(6));
          // 210
          off += fill_sym2(rule, off, Weight_(  8) / Weight_(60), Coord_(4) / Coord_(6), Coord_(1) / Coord_(6));
          // 111
          off += fill_sym3(rule, off, Weight_(  3) / Weight_(60), Coord_(3) / Coord_(6), Coord_(2) / Coord_(6), Coord_(1) / Coord_(6));
          break;

        case 4:
          // 400
          off += fill_sym2(rule, off, Weight_( 307) / Weight_(1440), Coord_(5) / Coord_(7), Coord_(1) / Coord_(7));
          // 211
          off += fill_sym2(rule, off, Weight_( -64) / Weight_(1440), Coord_(3) / Coord_(7), Coord_(2) / Coord_(7));
          // 022
          off += fill_sym2(rule, off, Weight_( 629) / Weight_(1440), Coord_(1) / Coord_(7), Coord_(3) / Coord_(7));
          // 310
          off += fill_sym3(rule, off, Weight_(-316) / Weight_(1440), Coord_(1) / Coord_(7), Coord_(2) / Coord_(7), Coord_(4) / Coord_(7));
          break;

        case 5:
          // TODO
          break;

        case 6:
          // TODO
          break;

        case 7:
          // TODO
          break;

        case 8:
          // 800
          off += fill_sym2(rule, off, Weight_( 1051445) / Weight_(7257600), Coord_(9) / Coord_(11), Coord_(1) / Coord_(11));
          // 611
          off += fill_sym2(rule, off, Weight_(-1818134) / Weight_(7257600), Coord_(7) / Coord_(11), Coord_(2) / Coord_(11));
          // 422
          off += fill_sym2(rule, off, Weight_(10685542) / Weight_(7257600), Coord_(5) / Coord_(11), Coord_(3) / Coord_(11));
          // 233
          off += fill_sym2(rule, off, Weight_(-6437608) / Weight_(7257600), Coord_(3) / Coord_(11), Coord_(4) / Coord_(11));
          // 044
          off += fill_sym2(rule, off, Weight_(12368047) / Weight_(7257600), Coord_(1) / Coord_(11), Coord_(5) / Coord_(11));
          // 710
          off += fill_sym3(rule, off, Weight_(-2366706) / Weight_(7257600), Coord_(8) / Coord_(11), Coord_(2) / Coord_(11), Coord_(1) / Coord_(11));
          // 620
          off += fill_sym3(rule, off, Weight_( 6493915) / Weight_(7257600), Coord_(7) / Coord_(11), Coord_(3) / Coord_(11), Coord_(1) / Coord_(11));
          // 530
          off += fill_sym3(rule, off, Weight_(-9986439) / Weight_(7257600), Coord_(6) / Coord_(11), Coord_(4) / Coord_(11), Coord_(1) / Coord_(11));
          // 521
          off += fill_sym3(rule, off, Weight_( 3757007) / Weight_(7257600), Coord_(6) / Coord_(11), Coord_(3) / Coord_(11), Coord_(2) / Coord_(11));
          // 431
          off += fill_sym3(rule, off, Weight_(  478257) / Weight_(7257600), Coord_(5) / Coord_(11), Coord_(4) / Coord_(11), Coord_(2) / Coord_(11));
          break;
        }
      }
    }; // class SilvesterOpenDriver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SILVESTER_DRIVER_HPP
