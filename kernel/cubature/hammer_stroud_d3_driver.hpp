#pragma once
#ifndef KERNEL_CUBATURE_HAMMER_STROUD_D3_DRIVER_HPP
#define KERNEL_CUBATURE_HAMMER_STROUD_D3_DRIVER_HPP 1

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

      class HammerStroudD3DriverBase :
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
          return "hammer-stroud-degree-3";
        }
      };
    } // namespace Intern
    /// \endcond

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
    class HammerStroudD3Driver;

    // Simplex specialisation
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class HammerStroudD3Driver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::HammerStroudD3DriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = dim_ + 2
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
        Weight_ B = - Weight_((dim_ + 1)*(dim_ + 1))/Weight_(4*dim_ + 8) * V;
        Weight_ C = Weight_((dim_ + 3)*(dim_ + 3))/Weight_(4*(dim_ + 1)*(dim_ + 2)) * V;

        // B-point
        rule.get_weight(0) = B;
        for(int j(0); j < dim_; ++j)
        {
          rule.get_coord(0,j) = Coord_(1)/Coord_(dim_ + 1);
        }

        // C-points
        Index count = 0;

        for(Index i(0); i <= Index(dim_); ++i)
        {
          ++count;

          // set weight
          rule.get_weight(count) = C;

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            if(Index(j) == i)
            {
              rule.get_coord(count,j) = Coord_(3)/Coord_(dim_ + 3);
            }
            else
            {
              rule.get_coord(count,j) = Coord_(1)/Coord_(dim_ + 3);
            }
          }
        }

      }
    }; // class HammerStroudD3Driver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_HAMMER_STROUD_D3_DRIVER_HPP