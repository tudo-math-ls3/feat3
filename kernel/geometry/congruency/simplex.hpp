#pragma once
#ifndef KERNEL_GEOMETRY_CONGRUENCY_SIMPLEX_HPP
#define KERNEL_GEOMETRY_CONGRUENCY_SIMPLEX_HPP 1

// includes, FEAST
#include <kernel/geometry/congruency/base.hpp>

namespace FEAST
{
  namespace Geometry
  {
    namespace Congruency
    {
      /**
       * \brief Congruency sampler implementation for Simplex<1> shape
         \verbatim
               0        1
             0---1    1---0
         \endverbatim
       * \author Peter Zajac
       */
      template<>
      class Sampler< Shape::Simplex<1> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          Index sv0 = src[0];
          if(sv0 == trg[0])
          {
            return 0;
          }
          else if(sv0 == trg[1])
          {
            return 1;
          }
          // invalid
          return -1;
        }
      }; // class Sampler<Simplex<1>>

      /**
       * \brief Congruency sampler implementation for Simplex<2> shape
         \verbatim
           0        1       2
              2        1       0
             / \      / \     / \
            0---1    2---0   1---2

           4        5       6
              1        2       0
             / \      / \     / \
            0---2    1---0   2---1
         \endverbatim
       *
       * \note The lack of the number \c 3 is intentional.
       *
       * \author Peter Zajac
       */
      template<>
      class Sampler< Shape::Simplex<2> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          Index sv0 = src[0];
          Index sv1 = src[1];
          if(sv0 == trg[0])
          {
            if(sv1 == trg[1])
            {
              return 0;
            }
            else if(sv1 == trg[2])
            {
              return 4;
            }
          }
          else if(sv0 == trg[1])
          {
            if(sv1 == trg[2])
            {
              return 1;
            }
            else if(sv1 == trg[0])
            {
              return 5;
            }
          }
          else if(sv0 == trg[2])
          {
            if(sv1 == trg[0])
            {
              return 2;
            }
            else if(sv1 == trg[1])
            {
              return 6;
            }
          }
          // invalid
          return -1;
        }
      }; // class Sampler<Simplex<2>>

      /// \cond internal
      /**
       * \brief Congruency vertex index mapping for Simplex<1> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class IndexMapping<Shape::Simplex<1>, 0>
      {
      public:
        static int map(int orient, int idx)
        {
          static const int indices[2][2] =
          {
            {0, 1},
            {1, 0}
          };

          return indices[orient][idx];
        }
      };

      /**
       * \brief Congruency vertex index mapping for Simplex<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class IndexMapping<Shape::Simplex<2>, 0>
      {
      public:
        static int map(int orient, int idx)
        {
          static const int indices[7][3] =
          {
            {0, 1, 2},
            {1, 2, 0},
            {2, 0, 1},
            {0, 0, 0}, // unused
            {0, 2, 1},
            {1, 0, 2},
            {2, 1, 0}
          };

          return indices[orient][idx];
        };
      };

      /**
       * \brief Congruency edge index mapping for Simplex<2> shape.
       *
       * \author Peter Zajac
       */

      template<>
      class IndexMapping<Shape::Simplex<2>, 1>
      {
      public:
        static int map(int orient, int idx)
        {
          static const int indices[7][3] =
          {
            {0, 1, 2},
            {1, 2, 0},
            {2, 0, 1},
            {0, 0, 0}, // unused
            {0, 2, 1},
            {1, 0, 2},
            {2, 1, 0}
          };

          return indices[orient][idx];
        };
      };

      /// \endcond
    } // namespace Congruency
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONGRUENCY_SIMPLEX_HPP
