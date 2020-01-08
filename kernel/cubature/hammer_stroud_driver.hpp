// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_HAMMER_STROUD_DRIVER_HPP
#define KERNEL_CUBATURE_HAMMER_STROUD_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/driver_base.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAT
{
  namespace Cubature
  {
    /**
     * \brief Hammer-Stroud-D2 driver class template
     *
     * This driver implements the Hammer-Stroud rule of degree two for simplices.
     * \see Stroud - Approximate Calculation Of Multiple Integrals,
     *      page 307, formula Tn:2-1
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Constantin Christof
     */
    template<typename Shape_>
    class HammerStroudD2Driver DOXY({});

    // Simplex specialisation
    template<int dim_>
    class HammerStroudD2Driver<Shape::Simplex<dim_> > :
      public DriverBase<Shape::Simplex<dim_> >
    {
    public:
      /// this rule is not variadic
      static constexpr bool variadic = false;
      static constexpr int num_points = dim_ + 1;

      ///Returns the name of the cubature rule.
      static String name()
      {
        return "hammer-stroud-degree-2";
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
      static void fill(Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_>& rule)
      {
        // auxiliary variables
        Coord_ r = (Coord_(dim_ + 2) - Math::sqrt(Coord_(dim_ + 2)))/
                   (Coord_(dim_ + 2) * Coord_(dim_ + 1));
        Coord_ s = (Coord_(dim_ + 2) + Coord_(dim_) * Math::sqrt(Coord_(dim_ + 2)))/
                   (Coord_(dim_ + 2) * Coord_(dim_ + 1));

        for(int i(0); i <= dim_; ++i)
        {
          // set weight
          rule.get_weight(i) = Weight_(1) / Weight_(MetaMath::Factorial<dim_ + 1>::value);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(i,j) = (i == j) ? s : r;
          }
        }
      }
    }; // class HammerStroudD2Driver<Simplex<...>,...>

    /**
     * \brief Hammer-Stroud-D3 driver class template
     *
     * This driver implements the Hammer-Stroud rule of degree three for simplices.
     * \see Stroud - Approximate Calculation Of Multiple Integrals,
     *      page 308, formula Tn:3-1
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Constantin Christof
     */
    template<typename Shape_>
    class HammerStroudD3Driver DOXY({});

    // Simplex specialisation
    template<int dim_>
    class HammerStroudD3Driver<Shape::Simplex<dim_> > :
      public DriverBase<Shape::Simplex<dim_> >
    {
    public:
      /// this rule is not variadic
      static constexpr bool variadic = false;
      static constexpr int num_points = dim_ + 2;

      ///Returns the name of the cubature rule.
      static String name()
      {
        return "hammer-stroud-degree-3";
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
      static void fill(Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_>& rule)
      {
        // auxiliary variables
        Weight_ V = Weight_(1) / Weight_(MetaMath::Factorial<dim_>::value);
        Weight_ B = - Weight_((dim_ + 1)*(dim_ + 1))/Weight_(4*dim_ + 8) * V;
        Weight_ C = Weight_((dim_ + 3)*(dim_ + 3))/Weight_(4*(dim_ + 1)*(dim_ + 2)) * V;

        // B-point
        rule.get_weight(0) = B;
        for(int j(0); j < dim_; ++j)
        {
          rule.get_coord(0,j) = Coord_(1) / Coord_(dim_ + 1);
        }

        // C-points
        int count = 0;

        for(int i(0); i <= dim_; ++i)
        {
          ++count;

          // set weight
          rule.get_weight(count) = C;

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            rule.get_coord(count,j) = Coord_((i == j) ? 3 : 1) / Coord_(dim_ + 3);
          }
        }

      }
    }; // class HammerStroudD3Driver<Simplex<...>,...>

    /**
     * \brief Hammer-Stroud-D5 driver class template
     *
     * This driver implements the Hammer-Stroud rule of degree five for tetrahedra.
     * \see Stroud - Approximate Calculation Of Multiple Integrals,
     *      page 315, formula T3:5-1
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Constantin Christof
     */
    template<typename Shape_>
    class HammerStroudD5Driver DOXY({});

    // Simplex specialisation
    template<>
    class HammerStroudD5Driver<Shape::Simplex<3> > :
      public DriverBase<Shape::Simplex<3> >
    {
    public:
      /// this rule is not variadic
      static constexpr bool variadic = false;
      static constexpr int dim = 3;
      static constexpr int num_points = 15;

      ///Returns the name of the cubature rule.
      static String name()
      {
        return "hammer-stroud-degree-5";
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
      static void fill(Rule<Shape::Simplex<3>, Weight_, Coord_, Point_>& rule)
      {
        // auxiliary variables
        Weight_ V = Weight_(1) / Weight_(MetaMath::Factorial<dim>::value);
        Weight_ A = Weight_(16)/Weight_(135) * V;
        Weight_ B1 = (Weight_(2665) + Weight_(14)*Math::sqrt(Weight_(15)))/Weight_(37800)*V;
        Weight_ B2 = (Weight_(2665) - Weight_(14)*Math::sqrt(Weight_(15)))/Weight_(37800)*V;
        Weight_ C = Weight_(20)/Weight_(378)*V;

        Coord_ s1 = (Coord_(7) - Math::sqrt(Coord_(15)))/Coord_(34);
        Coord_ s2 = (Coord_(7) + Math::sqrt(Coord_(15)))/Coord_(34);
        Coord_ t1 = (Coord_(13) + Coord_(3)*Math::sqrt(Coord_(15)))/Coord_(34);
        Coord_ t2 = (Coord_(13) - Coord_(3)*Math::sqrt(Coord_(15)))/Coord_(34);
        Coord_ u  = (Coord_(10) - Coord_(2)*Math::sqrt(Coord_(15)))/Coord_(40);
        Coord_ v  = (Coord_(10) + Coord_(2)*Math::sqrt(Coord_(15)))/Coord_(40);

        // counter
        int count = 0;

        // A-point
        rule.get_weight(count) = A;
        for(int j(0); j < dim; ++j)
        {
          rule.get_coord(count,j) = Coord_(1)/Coord_(4);
        }
        ++count;

        // B1-points
        for(int i(0); i <= dim; ++i)
        {
          // set weight
          rule.get_weight(count) = B1;

          // set point coords
          for(int j(0); j < dim; ++j)
          {
            rule.get_coord(count,j) = (i == j) ? t1 : s1;
          }
        ++count;
        }

        // B2-points
        for(int i(0); i <= dim; ++i)
        {

          // set weight
          rule.get_weight(count) = B2;

          // set point coords
          for(int j(0); j < dim; ++j)
          {
            rule.get_coord(count,j) = (i == j) ? t2 : s2;
          }
        ++count;
        }

        // C-points
        for(int i(0); i <= dim; ++i)
        {
          for(int j(0); j < i; ++j)
          {
            if(i != j)
            {
              // set weight
              rule.get_weight(count) = C;

              // set point coords
              for(int k(0); k < dim; ++k)
              {
                rule.get_coord(count,k) = (k == i) || (k == j) ? u : v;
              }
              ++count;
            }
          }
        }
      }
    }; // class HammerStroudD5Driver<Simplex<3>,...>
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_HAMMER_STROUD_DRIVER_HPP
