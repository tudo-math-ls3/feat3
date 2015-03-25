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
      /// this rule is variadic
      static constexpr bool variadic = true;
      static constexpr int min_points = 2;
      static constexpr int max_points = 8;

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
          // 111
          off += fill_sym1(rule, off, Weight_(-12) / Weight_(60), Coord_(2) / Coord_(6));
          // 300
          off += fill_sym2(rule, off, Weight_(  8) / Weight_(60), Coord_(4) / Coord_(6), Coord_(1) / Coord_(6));
          // 310
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
          // 500
          off += fill_sym2(rule, off, Weight_(  71) / Weight_(730), Coord_(6) / Coord_(8), Coord_(1) / Coord_(8));
          // 311
          off += fill_sym2(rule, off, Weight_(-167) / Weight_(730), Coord_(4) / Coord_(8), Coord_(2) / Coord_(8));
          // 122
          off += fill_sym2(rule, off, Weight_( 113) / Weight_(730), Coord_(2) / Coord_(8), Coord_(3) / Coord_(8));
          // 410
          off += fill_sym3(rule, off, Weight_( -13) / Weight_(730), Coord_(5) / Coord_(8), Coord_(2) / Coord_(8), Coord_(1) / Coord_(8));
          // 320
          off += fill_sym3(rule, off, Weight_(  57) / Weight_(730), Coord_(4) / Coord_(8), Coord_(3) / Coord_(8), Coord_(1) / Coord_(8));
          break;

        case 6:
          // 222
          off += fill_sym1(rule, off, Weight_( 3509) / Weight_(4480), Coord_(3) / Coord_(9));
          // 600
          off += fill_sym2(rule, off, Weight_(  767) / Weight_(4480), Coord_(7) / Coord_(9), Coord_(1) / Coord_(9));
          // 411
          off += fill_sym2(rule, off, Weight_(  387) / Weight_(4480), Coord_(4) / Coord_(9), Coord_(2) / Coord_(9));
          // 033
          off += fill_sym2(rule, off, Weight_(-3035) / Weight_(4480), Coord_(1) / Coord_(9), Coord_(4) / Coord_(9));
          // 510
          off += fill_sym3(rule, off, Weight_(-1257) / Weight_(4480), Coord_(6) / Coord_(9), Coord_(2) / Coord_(9), Coord_(1) / Coord_(9));
          // 420
          off += fill_sym3(rule, off, Weight_( 2901) / Weight_(4480), Coord_(5) / Coord_(9), Coord_(3) / Coord_(9), Coord_(1) / Coord_(9));
          // 321
          off += fill_sym3(rule, off, Weight_( -915) / Weight_(4480), Coord_(4) / Coord_(9), Coord_(3) / Coord_(9), Coord_(2) / Coord_(9));
          break;

        case 7:
          // 700
          off += fill_sym2(rule, off, Weight_(  898) / Weight_(9072), Coord_(8) / Coord_(10), Coord_(1) / Coord_(10));
          // 511
          off += fill_sym2(rule, off, Weight_(-2522) / Weight_(9072), Coord_(6) / Coord_(10), Coord_(2) / Coord_(10));
          // 133
          off += fill_sym2(rule, off, Weight_(-5726) / Weight_(9072), Coord_(2) / Coord_(10), Coord_(4) / Coord_(10));
          // 322
          off += fill_sym2(rule, off, Weight_( 1444) / Weight_(9072), Coord_(4) / Coord_(10), Coord_(3) / Coord_(10));
          // 610
          off += fill_sym3(rule, off, Weight_( -662) / Weight_(9072), Coord_(7) / Coord_(10), Coord_(2) / Coord_(10), Coord_(1) / Coord_(10));
          // 520
          off += fill_sym3(rule, off, Weight_( 1573) / Weight_(9072), Coord_(6) / Coord_(10), Coord_(3) / Coord_(10), Coord_(1) / Coord_(10));
          // 430
          off += fill_sym3(rule, off, Weight_( -191) / Weight_(9072), Coord_(5) / Coord_(10), Coord_(4) / Coord_(10), Coord_(1) / Coord_(10));
          // 421
          off += fill_sym3(rule, off, Weight_( 2989) / Weight_(9072), Coord_(5) / Coord_(10), Coord_(3) / Coord_(10), Coord_(2) / Coord_(10));
          break;

        case 8:
          // 800
          off += fill_sym2(rule, off, Weight_( 1051445) / Weight_(7257600), Coord_(9) / Coord_(11), Coord_(1) / Coord_(11));
          // 611
          off += fill_sym2(rule, off, Weight_( 1818134) / Weight_(7257600), Coord_(7) / Coord_(11), Coord_(2) / Coord_(11));
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
          // Note: The article referenced in the documentation of this class template claims that the following
          // weight should be positive rather than negative. However, in this case, the sum of the weights does
          // not match. The negative weight has been chosen in accordance to the implementation provided in
          // http://people.sc.fsu.edu/~jburkardt/f_src/nco_triangle/nco_triangle.html
          // 521
          off += fill_sym3(rule, off, Weight_(-3757007) / Weight_(7257600), Coord_(6) / Coord_(11), Coord_(3) / Coord_(11), Coord_(2) / Coord_(11));
          // 431
          off += fill_sym3(rule, off, Weight_(  478257) / Weight_(7257600), Coord_(5) / Coord_(11), Coord_(4) / Coord_(11), Coord_(2) / Coord_(11));
          break;
        }
      }
    }; // class SilvesterOpenDriver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SILVESTER_DRIVER_HPP
