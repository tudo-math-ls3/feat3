#pragma once
#ifndef KERNEL_CUBATURE_LAUFFER_D2_DRIVER_HPP
#define KERNEL_CUBATURE_LAUFFER_D2_DRIVER_HPP 1

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

      class LaufferD2DriverBase :
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
          return "lauffer-degree-2";
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Lauffer-D2 driver class template
     *
     * This driver implements the Lauffer rule of degree two for simplices.
     * \see Stroud - Approximate Calculation Of Multiple Integrals,
     *      page 307, formula Tn:2-2
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
    class LaufferD2Driver;

    // Simplex specialisation
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class LaufferD2Driver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::LaufferD2DriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = (dim_ + 1)*(dim_ + 2)/2
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
        Weight_ B = Weight_(2 - dim_)/Weight_((dim_ + 1)*(dim_ + 2)) * V;
        Weight_ C = Weight_(4)/Weight_((dim_ + 1)*(dim_ + 2)) * V;

        // B-points
        for(Index i(0); i <= Index(dim_); ++i)
        {
          // set weight
          rule.get_weight(i) = B;

          // set point coords
          for(int j(0); j < dim_; ++j)
          {
            if(Index(j) == i)
            {
              rule.get_coord(i,j) = 1;
            }
            else
            {
              rule.get_coord(i,j) = 0;
            }
          }
        }

        // counter
        Index count = Index(dim_);

        // C-points
        for(Index i(1); i <= Index(dim_); ++i)
        {
          for(Index k(0); k < i; ++k)
          {
            ++count;
            rule.get_weight(count) = C;

            // set point coords
            for(int j(0); j < dim_; ++j)
            {
              if(Index(j) == k || Index(j) == i)
              {
                rule.get_coord(count,j) = Coord_(1)/Coord_(2);
              }
              else
              {
                rule.get_coord(count,j) = 0;
              }
            } // j-loop
          } // k-loop
        } // i-loop

      }
    }; // class LaufferD2Driver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_LAUFFER_D2_DRIVER_HPP