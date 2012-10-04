#pragma once
#ifndef KERNEL_CUBATURE_HAMMER_STROUD_D2_DRIVER_HPP
#define KERNEL_CUBATURE_HAMMER_STROUD_D2_DRIVER_HPP 1

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
      class HammerStroudD2DriverBase :
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
          return "hammer-stroud-degree-2";
        }
      };
    } // namespace Intern
    /// \endcond

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
    class HammerStroudD2Driver;

    // Simplex specialisation
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class HammerStroudD2Driver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::HammerStroudD2DriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
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
      static void fill(RuleType& rule)
      {
        // auxiliary variables
        Coord_ r = (Coord_(dim_ + 2) - std::sqrt(Coord_(dim_ + 2)))/
                   (Coord_(dim_ + 2) * Coord_(dim_ + 1));
        Coord_ s = (Coord_(dim_ + 2) + Coord_(dim_) * std::sqrt(Coord_(dim_ + 2)))/
                   (Coord_(dim_ + 2) * Coord_(dim_ + 1));

        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = Weight_(1) / Weight_(Factorial<dim_ + 1>::value);

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            if(Index(j) == i)
            {
              rule.get_coord(i,j) = s;
            }
            else
            {
              rule.get_coord(i,j) = r;
            }
          }
        }
      }
    }; // class HammerStroudD2Driver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_HAMMER_STROUD_D2_DRIVER_HPP