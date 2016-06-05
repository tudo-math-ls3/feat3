#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_FACE_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_FACE_INDEX_MAPPING_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Face index mapping class template
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int cell_dim_,
        int face_dim_>
#ifndef DOXYGEN
      class FaceIndexMapping;
#else
      class FaceIndexMapping
      {
      public:
        static int map(int cell, int face);
      };
#endif // DOXYGEN

      template<>
      class FaceIndexMapping<Shape::Simplex<2>, 1, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 3), "invalid cell index");
          ASSERTM((0 <= face) && (face < 2), "invalid face index");
          static const int indices[3][2] =
          {
            {1, 2},
            {2, 0},
            {0, 1}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Simplex<2>, 1, 0>

      template<>
      class FaceIndexMapping<Shape::Simplex<3>, 1, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 6), "invalid cell index");
          ASSERTM((0 <= face) && (face < 2), "invalid face index");
          static const int indices[6][2] =
          {
            {0, 1},
            {0, 2},
            {0, 3},
            {1, 2},
            {1, 3},
            {2, 3}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Simplex<3>, 1, 0>

      template<>
      class FaceIndexMapping<Shape::Simplex<3>, 2, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 4), "invalid cell index");
          ASSERTM((0 <= face) && (face < 3), "invalid face index");
          static const int indices[4][3] =
          {
            {1, 2, 3},
            {0, 2, 3},
            {0, 1, 3},
            {0, 1, 2}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Simplex<3>, 2, 0>

      template<>
      class FaceIndexMapping<Shape::Simplex<3>, 2, 1>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 4), "invalid cell index");
          ASSERTM((0 <= face) && (face < 3), "invalid face index");
          static const int indices[4][3] =
          {
            {5, 4, 3},
            {5, 2, 1},
            {4, 2, 0},
            {3, 1, 0}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Simplex<3>, 2, 0>

      template<>
      class FaceIndexMapping<Shape::Hypercube<2>, 1, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 4), "invalid cell index");
          ASSERTM((0 <= face) && (face < 2), "invalid face index");
          static const int indices[4][2] =
          {
            {0, 1},
            {2, 3},
            {0, 2},
            {1, 3}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Hypercube<2>, 1, 0>

      template<>
      class FaceIndexMapping<Shape::Hypercube<3>, 1, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 12), "invalid cell index");
          ASSERTM((0 <= face) && (face < 2), "invalid face index");
          static const int indices[12][2] =
          {
            {0, 1},
            {2, 3},
            {4, 5},
            {6, 7},
            {0, 2},
            {1, 3},
            {4, 6},
            {5, 7},
            {0, 4},
            {1, 5},
            {2, 6},
            {3, 7}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Shape::Hypercube<3>, 1, 0>

      template<>
      class FaceIndexMapping<Shape::Hypercube<3>, 2, 0>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 6), "invalid cell index");
          ASSERTM((0 <= face) && (face < 4), "invalid face index");
          static const int indices[6][4] =
          {
            {0, 1, 2, 3},
            {4, 5, 6, 7},
            {0, 1, 4, 5},
            {2, 3, 6, 7},
            {0, 2, 4, 6},
            {1, 3, 5, 7}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Shape::Hypercube<3>, 2, 0>

      template<>
      class FaceIndexMapping<Shape::Hypercube<3>, 2, 1>
      {
      public:
        static int map(int cell, int face)
        {
          ASSERTM((0 <= cell) && (cell < 6), "invalid cell index");
          ASSERTM((0 <= face) && (face < 4), "invalid face index");
          static const int indices[6][4] =
          {
            {0, 1,  4,  5},
            {2, 3,  6,  7},
            {0, 2,  8,  9},
            {1, 3, 10, 11},
            {4, 6,  8, 10},
            {5, 7,  9, 11}
          };
          return indices[cell][face];
        }
      }; // FaceIndexMapping<Shape::Hypercube<3>, 2, 1>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_FACE_INDEX_MAPPING_HPP
