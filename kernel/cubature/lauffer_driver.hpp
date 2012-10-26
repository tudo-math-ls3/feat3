#pragma once
#ifndef KERNEL_CUBATURE_LAUFFER_DRIVER_HPP
#define KERNEL_CUBATURE_LAUFFER_DRIVER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_base.hpp>

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
    class LaufferD2Driver DOXY({});

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
      static void fill(Rule<Shape::Simplex<dim_>, Weight_, Coord_, Point_>& rule)
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
    class LaufferD4Driver DOXY({});

    // Tetrahedron specialisation
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class LaufferD4Driver<Shape::Simplex<3>, Weight_, Coord_, Point_> :
      public Intern::LaufferD4DriverBase
    {
    public:
      enum
      {
        dim = 3,
        num_points = (Factorial<dim + 4>::value)/(24*Factorial<dim>::value)
      };

      /**
       * \brief Fills the cubature rule structure.
       *
       * \param[in,out] rule
       * The cubature rule to be filled.
       */
      static void fill(Rule<Shape::Simplex<3>, Weight_, Coord_, Point_>& rule)
      {
        // auxiliary variables
        Weight_ V = Weight_(1) / Weight_(Factorial<dim>::value);
        Weight_ B1 =
          Weight_(-3*dim*dim*dim + 17*dim*dim - 58*dim + 72)/
          Weight_(3*(dim + 1)*(dim + 2)*(dim + 3)*(dim + 4)) * V;
        Weight_ B2 =
          Weight_(16*(dim*dim - 5*dim + 12))/
          Weight_(3*(dim + 1)*(dim + 2)*(dim + 3)*(dim + 4)) * V;
        Weight_ B3 =
          Weight_(4*(dim*dim - 9*dim + 12))/
          Weight_((dim + 1)*(dim + 2)*(dim + 3)*(dim + 4)) * V;
        Weight_ B4 =
          Weight_(64*(4 - dim))/
          Weight_(2*(dim + 1)*(dim + 2)*(dim + 3)*(dim + 4)) * V;
        Weight_ B5 =
          Weight_(256)/
          Weight_((dim + 1)*(dim + 2)*(dim + 3)*(dim + 4)) * V;

        Index count = 0;

        // B1-points
        for(Index i(0); i <= Index(dim); ++i)
        {
          // set weight
          rule.get_weight(count) = B1;

          // set point coords
          for(int j(0); j < dim; ++j)
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
        for(Index i(0); i <= Index(dim); ++i)
        {
          for(Index j(0); j <= Index(dim); ++j)
          {
            if(i != j)
            {
              // set weight
              rule.get_weight(count) = B2;

              // set point coords
              for(int k(0); k < dim; ++k)
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
        for(Index i(1); i <= Index(dim); ++i)
        {
          for(Index j(0); j < i; ++j)
          {
            rule.get_weight(count) = B3;

            // set point coords
            for(int k(0); k < dim; ++k)
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
        for(Index i(0); i <= Index(dim); ++i)
        {
          for(Index j(0); j <= Index(dim); ++j)
          {
            for(Index k(0); k < j; ++k)
            {
              if(i != j && i != k && j != k)
              {
                rule.get_weight(count) = B4;

                // set point coords
                for(int l(0); l < dim; ++l)
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
        for(Index i(1); i <= Index(dim); ++i)
        {
          for(Index j(0); j < i; ++j)
          {
            for(Index k(0); k < j; ++k)
            {
              for(Index l(0); l < k; ++l)
              {
                rule.get_weight(count) = B5;

                // set point coords
                for(int m(0); m < dim; ++m)
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
    }; // class LaufferD4Driver<Simplex<3>,...>

  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_LAUFFER_DRIVER_HPP