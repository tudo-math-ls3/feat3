#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_CONGRUENCY_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_CONGRUENCY_MAPPING_HPP 1

// includes, FEAT
#include <kernel/geometry/intern/congruency_sampler.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Congruency index mapping class template.
       *
       * \todo detailed description
       *
       * \tparam Shape_
       * The shape tag class for the index mapping.
       *
       * \tparam face_dim_
       * The dimension of the faces for the index mapping.
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int face_dim_>
#ifndef DOXYGEN
      class CongruencyMapping;
#else
      class CongruencyMapping
      {
      public:
        /**
         * \brief Maps a local index based on an orientation code.
         *
         * \param[in] orient
         * An orientation code as returned by the CongruencySampler<Shape_>::compare() function.
         *
         * \param[in] idx
         * The local face index that is to be mapped.
         *
         * \returns
         * The mapped local face index.
         */
        static int map(int orient, int idx);
      };
#endif // DOXYGEN

      /**
       * \brief Congruency vertex index mapping for Simplex<1> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencyMapping<Shape::Simplex<1>, 0>
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
      class CongruencyMapping<Shape::Simplex<2>, 0>
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
      class CongruencyMapping<Shape::Simplex<2>, 1>
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
       * \brief Congruency vertex index mapping for Hypercube<1> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencyMapping<Shape::Hypercube<1>, 0>
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
       * \brief Congruency vertex index mapping for Hypercube<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencyMapping<Shape::Hypercube<2>, 0>
      {
      public:
        static int map(int orient, int idx)
        {
          static const int indices[8][4] =
          {
            {0, 1, 2, 3},
            {1, 3, 0, 2},
            {2, 0, 3, 1},
            {3, 2, 1, 0},
            {0, 2, 1, 3},
            {1, 0, 3, 2},
            {2, 3, 0, 1},
            {3, 1, 2, 0}
          };

          return indices[orient][idx];
        }
      };

      /**
       * \brief Congruency edge index mapping for Hypercube<2> shape.
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencyMapping<Shape::Hypercube<2>, 1>
      {
      public:
        static int map(int orient, int idx)
        {
          static const int indices[8][4] =
          {
            {0, 1, 2, 3},
            {3, 2, 0, 1},
            {2, 3, 1, 0},
            {1, 0, 3, 2},
            {2, 3, 0, 1},
            {0, 1, 3, 2},
            {1, 0, 2, 3},
            {3, 2, 1, 0}
          };

          return indices[orient][idx];
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_CONGRUENCY_MAPPING_HPP
