// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_CODING_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_CODING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<int shape_dim_, int cell_dim_>
      struct StructIndexCoding;

      // 1D vertices
      template<>
      struct StructIndexCoding<1, 0>
      {
        // n[0] = number of slices in X
        // v[0] = x-index of vertex <= n[0]
        static Index encode(const Index v[], const Index /*n*/[])
        {
          return v[0];
        }

        static void decode(Index v[], const Index i, const Index /*n*/[])
        {
          v[0] = i;
        }
      };

      // 2D vertices
      template<>
      struct StructIndexCoding<2, 0>
      {
        // n[0] = number of slices in X
        // n[1] = number of slices in Y
        // v[0] = x-index of vertex <= n[0]
        // v[1] = y-index of vertex <= n[1]
        static Index encode(const Index v[], const Index n[])
        {
          return v[0] + (n[0] + Index(1)) * v[1];
        }

        static void decode(Index v[], const Index i, const Index n[])
        {
          v[0] = i % (n[0] + Index(1));
          v[1] = i / (n[0] + Index(1));
        }
      };

      // 3D vertices
      template<>
      struct StructIndexCoding<3, 0>
      {
        // n[0] = number of slices in X
        // n[1] = number of slices in Y
        // n[2] = number of slices in Z
        // v[0] = x-index of vertex <= n[0]
        // v[1] = y-index of vertex <= n[1]
        // v[2] = z-index of vertex <= n[2]
        static Index encode(const Index v[], const Index n[])
        {
          return v[0] + (n[0] + Index(1)) * (v[1] + (n[1] + Index(1)) * v[2]);
        }

        static void decode(Index v[], const Index i, const Index n[])
        {
          v[0] = i % ((n[0] + Index(1)) * (n[1] + Index(1)));
          v[1] = (i / (n[0] + Index(1))) % (n[1] + Index(1));
          v[2] = i / ((n[0] + Index(1)) * (n[1] + Index(1)));
        }
      };

      // 1D edges
      template<>
      struct StructIndexCoding<1, 1>
      {
        // n[0] = number of slices in X
        // v[1] = x-index of edge < n[0]
        static Index encode(const Index v[], const Index /*n*/[])
        {
          return v[0];
        }

        static void decode(Index v[], const Index i, const Index /*n*/[])
        {
          v[0] = i;
        }
      };

      // 2D edges
      template<>
      struct StructIndexCoding<2, 1>
      {
        // n[0] = number of slices in X
        // n[1] = number of slices in Y
        // v[0] = x-index of edge < n[0] (+1)
        // v[1] = y-index of edge < n[1] (+1)
        // v[2] = {0 = horz edge, 1 = vert edge}
        static Index encode(const Index v[], const Index n[])
        {
          if(v[2] == 0) // horizontal edge
            return v[0] + n[0] * v[1];
          else // vertical edge
            return v[1] + n[1] * v[0] + n[0]*(n[1]+Index(1));
        }

        static void decode(Index v[], const Index i, const Index n[])
        {
          if(i < n[0] * (n[1] + Index(1)))
          {
            // horizontal edge
            v[0] = i % n[0];
            v[1] = i / n[0];
            v[2] = 0;
          }
          else
          {
            // vertical edge
            const Index j = (i - n[0] * (n[1] + Index(1)));
            v[0] = j / n[1];
            v[1] = j % n[1];
            v[2] = 1;
          }
        }
      };

      /// \todo 3D edges

      // 2D quads
      template<>
      struct StructIndexCoding<2, 2>
      {
        // n[0] = number of slices in X
        // n[1] = number of slices in Y
        // v[0] = x-index of quad < n[0]
        // v[1] = y-index of quad < n[1]
        static Index encode(const Index v[], const Index n[])
        {
          return v[0] + n[0] * v[1];
        }

        static void decode(Index v[], const Index i, const Index n[])
        {
          v[0] = i % n[0];
          v[1] = i / n[0];
        }
      };

      /// \todo 3D quads

      // 3D hexas
      template<>
      struct StructIndexCoding<3, 3>
      {
        // n[0] = number of slices in X
        // n[1] = number of slices in Y
        // n[2] = number of slices in Z
        // v[0] = x-index of hexa < n[0]
        // v[1] = y-index of hexa < n[1]
        // v[2] = z-index of hexa < n[2]
        static Index encode(const Index v[], const Index n[])
        {
          return v[0] + n[0] * (v[1] + n[1] * v[2]);
        }

        static void decode(Index v[], const Index i, const Index n[])
        {
          v[0] = i % (n[0] * n[1]);
          v[1] = (i / n[0]) % (n[1]);
          v[2] = i / (n[0] * n[1]);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_CODING_HPP
