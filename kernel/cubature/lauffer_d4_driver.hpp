#pragma once
#ifndef KERNEL_CUBATURE_LAUFFER_D4_DRIVER_HPP
#define KERNEL_CUBATURE_LAUFFER_D4_DRIVER_HPP 1

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
      class LaufferD4DriverBase :
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
          return "lauffer-degree-4";
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Lauffer-D4 driver class template
     *
     * This driver implements the Lauffer rule of degree four for tetrahedra.
     * \see Stroud - Approximate Calculation Of Multiple Integrals,
     *      page 311, formula Tn:4-1
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
    class LaufferD4Driver;

    // Simplex specialisation
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class LaufferD4Driver<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public Intern::LaufferD4DriverBase
    {
    public:
      typedef Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_> RuleType;
      enum
      {
        num_points = (Factorial<dim_ + 4>::value)/(24*Factorial<dim_>::value)
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
        Weight_ B1 =
          Weight_(-3*dim_*dim_*dim_ + 17*dim_*dim_ - 58*dim_ + 72)/
          Weight_(3*(dim_ + 1)*(dim_ + 2)*(dim_ + 3)*(dim_ + 4)) * V;
        Weight_ B2 =
          Weight_(16*(dim_*dim_ - 5*dim_ + 12))/
          Weight_(3*(dim_ + 1)*(dim_ + 2)*(dim_ + 3)*(dim_ + 4)) * V;
        Weight_ B3 =
          Weight_(4*(dim_*dim_ - 9*dim_ + 12))/
          Weight_((dim_ + 1)*(dim_ + 2)*(dim_ + 3)*(dim_ + 4)) * V;
        Weight_ B4 =
          Weight_(64*(4 - dim_))/
          Weight_(2*(dim_ + 1)*(dim_ + 2)*(dim_ + 3)*(dim_ + 4)) * V;
        Weight_ B5 =
          Weight_(256)/
          Weight_((dim_ + 1)*(dim_ + 2)*(dim_ + 3)*(dim_ + 4)) * V;

        Index count = 0;

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
              rule.get_coord(count,j) = 1;
            }
            else
            {
              rule.get_coord(count,j) = 0;
            }
          }
          ++count;
        }

        // B2-points
        for(Index i(0); i <= Index(dim_); ++i)
        {
          for(Index j(0); j <= Index(dim_); ++j)
          {
            if(i != j)
            {
              // set weight
              rule.get_weight(count) = B2;

              // set point coords
              for(int k(0); k < dim_; ++k)
              {
                if(Index(k) == i)
                {
                  rule.get_coord(count,k) = Coord_(1)/Coord_(4);
                }
                else if (Index(k) == j)
                {
                  rule.get_coord(count,k) = Coord_(3)/Coord_(4);
                }
                else
                {
                  rule.get_coord(count,k) = Coord_(0);
                }
              }
              ++count;
            }
          }
        }

        // B3-points
        for(Index i(1); i <= Index(dim_); ++i)
        {
          for(Index j(0); j < i; ++j)
          {
            rule.get_weight(count) = B3;

            // set point coords
            for(int k(0); k < dim_; ++k)
            {
              if(Index(k) == j || Index(k) == i)
              {
                rule.get_coord(count,k) = Coord_(1)/Coord_(2);
              }
              else
              {
                rule.get_coord(count,k) = 0;
              }
            } // k-loop
            ++count;
          } // j-loop
        } // i-loop

        // B4-points
        for(Index i(0); i <= Index(dim_); ++i)
        {
          for(Index j(0); j <= Index(dim_); ++j)
          {
            for(Index k(0); k < j; ++k)
            {
              if(i != j && i != k && j != k)
              {
                rule.get_weight(count) = B4;

                // set point coords
                for(int l(0); l < dim_; ++l)
                {
                  if(Index(l) == j || Index(l) == k)
                  {
                    rule.get_coord(count,l) = Coord_(1)/Coord_(4);
                  }
                  else if(Index(l) == i)
                  {
                    rule.get_coord(count,l) = Coord_(1)/Coord_(2);
                  }
                  else
                  {
                    rule.get_coord(count,l) = Coord_(0);
                  }
                } // l-loop
                ++count;
              } // if
            } // k-loop
          } // j-loop
        } //i-loop

        // B5-points
        for(Index i(1); i <= Index(dim_); ++i)
        {
          for(Index j(0); j < i; ++j)
          {
            for(Index k(0); k < j; ++k)
            {
              for(Index l(0); l < k; ++l)
              {
                rule.get_weight(count) = B5;

                // set point coords
                for(int m(0); m < dim_; ++m)
                {
                  if(Index(m) == i || Index(m) == j || Index(m) == k || Index(m) == l)
                  {
                    rule.get_coord(count,m) = Coord_(1)/Coord_(4);
                  }
                  else
                  {
                    rule.get_coord(count,m) = Coord_(0);
                  }
                } // m-loop
                ++count;
              } // l-loop
            } // k-loop
          } // j-loop
        } // i-loop
      }

    }; // class LaufferD4Driver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_LAUFFER_D4_DRIVER_HPP