#pragma once
#ifndef KERNEL_CUBATURE_HAMMER_STROUD_D5_DRIVER_HPP
#define KERNEL_CUBATURE_HAMMER_STROUD_D5_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_base.hpp>

// includes, STL
#include <cmath>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {

      class HammerStroudD5DriverBase :
        public DriverBase
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
          return "hammer-stroud-degree-5";
        }
      };
    } // namespace Intern
    /// \endcond

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
     * \tparam Weight_
     * The data type of the cubature weights.
     *
     * \tparam Coord_
     * The data type of the cubature point coordinates.
     *
     * \tparam Point_
     *
     * \author Constantin Christof
     */
    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class HammerStroudD5Driver;

    // Simplex specialisation
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class HammerStroudD5Driver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::HammerStroudD5DriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = 15
      };

      /**
       * \brief Fills the cubature rule structure.
       *
       * \param[in,out] rule
       * The cubature rule to be filled.
       */
      static void fill(RuleType& rule)
      {
        // auxiliary variables
        Weight_ V = Weight_(1) / Weight_(Factorial<dim_>::value);
        Weight_ A = Weight_(16)/Weight_(135) * V;
        Weight_ B1 = (Weight_(2665) + Weight_(14)*std::sqrt(Weight_(15)))/Weight_(37800)*V;
        Weight_ B2 = (Weight_(2665) - Weight_(14)*std::sqrt(Weight_(15)))/Weight_(37800)*V;
        Weight_ C = Weight_(20)/Weight_(378)*V;

        Coord_ s1 = (Coord_(7) - std::sqrt(Coord_(15)))/Coord_(34);
        Coord_ s2 = (Coord_(7) + std::sqrt(Coord_(15)))/Coord_(34);
        Coord_ t1 = (Coord_(13) + Coord_(3)*std::sqrt(Coord_(15)))/Coord_(34);
        Coord_ t2 = (Coord_(13) - Coord_(3)*std::sqrt(Coord_(15)))/Coord_(34);
        Coord_ u  = (Coord_(10) - Coord_(2)*std::sqrt(Coord_(15)))/Coord_(40);
        Coord_ v  = (Coord_(10) + Coord_(2)*std::sqrt(Coord_(15)))/Coord_(40);

        // counter
        Index count = 0;

        // A-point
        rule.get_weight(count) = A;
        for(int j(0); j < dim_; ++j)
        {
          rule.get_coord(count,j) = Coord_(1)/Coord_(4);
        }
        ++count;

        // B1-points
        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(count) = B1;

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            if(Index(j) == i)
            {
              rule.get_coord(count,j) = t1;
            }
            else
            {
              rule.get_coord(count,j) = s1;
            }
          }
        ++count;
        }

        // B2-points
        for(Index i(0); i <= Index(dim_); ++i)
        {

          // set weight
          rule.get_weight(count) = B2;

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            if(Index(j) == i)
            {
              rule.get_coord(count,j) = t2;
            }
            else
            {
              rule.get_coord(count,j) = s2;
            }
          }
        ++count;
        }

        // C-points
        for(Index i(0); i <= Index(dim_); ++i)
        {
          for(Index j(0); j < i; ++j)
          {
            if(i != j)
            {
              // set weight
              rule.get_weight(count) = C;

              // set point coords
              for(int k(0); k < dim_; ++k)
              {
                if(Index(k) == i || Index(k) == j)
                {
                  rule.get_coord(count,k) = u;
                }
                else
                {
                  rule.get_coord(count,k) = v;
                }
              }
              ++count;
            }
          }
        }
      }
    }; // class HammerStroudD5Driver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_HAMMER_STROUD_D5_DRIVER_HPP