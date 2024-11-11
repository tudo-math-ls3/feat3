// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/intern/congruency_sampler.hpp>

// includes, system
#include <algorithm>

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

        /**
         * \brief Flips the indices of the shape to invert the orientation
         *
         * \param[inout] idx
         * The index tuple which is to be flipped to invert its orientation
         */
        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx);
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

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // swap vertex 0 and 1
          std::swap(idx[0], idx[1]);
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
        }

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // 2          1
          // |'\.   --> |'\.
          // 0---1      0---2
          std::swap(idx[1], idx[2]);
        }
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
        }

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // 2          1
          // |'\.   --> |'\.
          // 0---1      0---2
          //
          // X          X
          // 2'1.   --> 0'1.
          // X-0-X      X-2-X
          //
          // swap edge 0 and 2
          std::swap(idx[0], idx[2]);
        }
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

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // swap vertex 0 and 1
          std::swap(idx[0], idx[1]);
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

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // 2---3      0---1
          // |   |  --> |   |
          // 0---1      2---3
          std::swap(idx[0], idx[2]);
          std::swap(idx[1], idx[3]);
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

        template<typename IdxTuple_>
        static void flip(IdxTuple_& idx)
        {
          // 2---3      0---1
          // |   |  --> |   |
          // 0---1      2---3
          //
          // X-1-X      X-0-X
          // 2   3  --> 2   3
          // X-0-X      X-1-X
          //
          // swap edge 0 and 1
          std::swap(idx[0], idx[1]);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
