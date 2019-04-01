// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SYMMETRIC_SIMPLEX_DRIVER_HPP
#define KERNEL_CUBATURE_SYMMETRIC_SIMPLEX_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/driver_base.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAT
{
  namespace Cubature
  {
    /**
     * \brief Symmetric Simplex Driver helper class.
     *
     * This class acts as a helper for symmetric simplex cubature drivers.
     *
     * \todo implement for Simplex<3>
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class SymmetricSimplexDriver DOXY({});

    // Simplex specialisation
    template<>
    class SymmetricSimplexDriver<Shape::Simplex<2> > :
      public DriverBase<Shape::Simplex<2> >
    {
    protected:
      /// Adds the one permutation of the barycentric point (x0,x0,x0) to the rule.
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static int fill_sym1(
        Rule<Shape::Simplex<2>, Weight_, Coord_, Point_>& rule,
        int off,
        Weight_ w,
        Coord_ x0)
      {
        rule.get_weight(off) = w;
        rule.get_coord(off, 0) = x0; // 000
        rule.get_coord(off, 1) = x0;
        return 1;
      }

      /// Adds the three permutations of the barycentric point (x0,x1,x1) to the rule.
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static int fill_sym2(
        Rule<Shape::Simplex<2>, Weight_, Coord_, Point_>& rule,
        int off,
        Weight_ w,
        Coord_ x0,
        Coord_ x1)
      {
        rule.get_weight(off) = w;
        rule.get_coord(off, 0) = x0; // 011
        rule.get_coord(off, 1) = x1;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x1; // 101
        rule.get_coord(off, 1) = x0;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x1; // 110
        rule.get_coord(off, 1) = x1;
        return 3;
      }

      /// Adds the six permutations of the barycentric point (x0,x1,x2) to the rule.
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static int fill_sym3(
        Rule<Shape::Simplex<2>, Weight_, Coord_, Point_>& rule,
        int off,
        Weight_ w,
        Coord_ x0,
        Coord_ x1,
        Coord_ x2)
      {
        rule.get_weight(off) = w;
        rule.get_coord(off, 0) = x0; // 012
        rule.get_coord(off, 1) = x1;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x1; // 120
        rule.get_coord(off, 1) = x2;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x2; // 201
        rule.get_coord(off, 1) = x0;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x1; // 102
        rule.get_coord(off, 1) = x0;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x2; // 210
        rule.get_coord(off, 1) = x1;
        rule.get_weight(++off) = w;
        rule.get_coord(off, 0) = x0; // 021
        rule.get_coord(off, 1) = x2;
        return 6;
      }
    }; // class SymmetricSimplexDriver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SYMMETRIC_SIMPLEX_DRIVER_HPP
