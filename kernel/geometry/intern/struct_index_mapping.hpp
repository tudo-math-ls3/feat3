// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        int shape_dim_,
        int cell_dim_,
        int face_dim_>
      struct StructIndexMapping;


      // one-dimensional

      /// vert-at-edge for edge mesh
      template<>
      struct StructIndexMapping<1, 1, 0>
      {
        static Index compute(Index i, Index j, const Index /*num_slices*/[])
        {
          return i + j;
        }
      }; //<1, 1, 0>


      // quad mesh

      /// vert-at-edge for quad mesh
      template<>
      struct StructIndexMapping<2, 1, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          if(i < num_slices[0] * (num_slices[1] + 1))
          {
            return i + j + i / num_slices[0];
          }
          else
          {
            return ((i - num_slices[0] * (num_slices[1] + 1)) / num_slices[1])
                   + ((i - num_slices[0] * (num_slices[1] + 1)) % num_slices[1])
                   * (num_slices[0] + 1) + j * (num_slices[0] + 1);
          }
        }
      }; //<2, 1, 0>

      /// vert-at-quad for quad mesh
      template<>
      struct StructIndexMapping<2, 2, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          return ((i / num_slices[0]) + ((j >> 1) & 0x1)) * (num_slices[0] + 1) + (i % num_slices[0]) + (j & 0x1);
        }
      }; //<2, 2, 0>

      /// edge-at-quad for quad mesh
      template<>
      struct StructIndexMapping<2, 2, 1>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          return ((3 - j) / 2) * (i + j * num_slices[0]) + (j / 2) * (num_slices[0] * (num_slices[1] + 1)
                  + ((i % num_slices[0]) + (j / 3)) * num_slices[1] + (i / num_slices[0]));
        }
      }; //<2, 2, 1>


      // hexahedron mesh

      /// vert-at-edge for hexa mesh
      template<>
      struct StructIndexMapping<3, 1, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // auxiliary index
          Index pos;

          // compute auxiliary variable

          //z-direction
          if( i >= num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
            (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
          {
            pos = i - num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) - (num_slices[2] + 1)
                  * (num_slices[1] + 1) * num_slices[0];
            return (pos / num_slices[2]) + (pos % num_slices[2] + j) * (num_slices[0] + 1)*(num_slices[1] + 1);
          }
          //y-direction
          else if( i >= (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0])
          {
            pos = i - (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0];
            Index x = (pos / num_slices[1]) % (num_slices[0] + 1);
            Index y = pos % num_slices[1];
            Index z = (pos / num_slices[1]) / (num_slices[0] + 1);
            return z * (num_slices[0] + 1) * (num_slices[1] + 1) + (y + j) * (num_slices[0] + 1) + x;
          }
          // x-direction
          else
          {
          return i + j + i / num_slices[0];
          }
        }
      }; //<3, 1, 0>

      /// vert-at-quad for hexahedron mesh
      template<>
      struct StructIndexMapping<3, 2, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          Index pos;
          if( i < num_slices[0] * num_slices[1] * (num_slices[2] + 1))
          {
            return i + (i / num_slices[0]) + (i / (num_slices[0] * num_slices[1])) * (num_slices[0] + 1)
                   + (j % 2) + (j / 2) * (num_slices[0] + 1);
          }
          else if( i >= num_slices[0] * num_slices[1] * (num_slices[2] + 1)
                   && i < num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1)))
          {
            pos = i - num_slices[0] * num_slices[1] * (num_slices[2] + 1);
            return pos % num_slices[0] + (pos / (num_slices[0] * num_slices[2])) * (num_slices[0] + 1)
                   + (pos % (num_slices[0] * num_slices[2])) / num_slices[0]
                   *(num_slices[0] + 1)*(num_slices[1] + 1) + (j % 2) + (j / 2) * (num_slices[0] + 1) * (num_slices[1] + 1);
          }
          else
          {
            pos = i - num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1));
            return pos / (num_slices[1] * num_slices[2]) + (pos % num_slices[1]) * (num_slices[0] + 1) + ((pos / num_slices[1]) % num_slices[2])
                   * (num_slices[0] + 1) * (num_slices[1] + 1) + (j % 2) * (num_slices[0] + 1) + (j / 2)
                   * (num_slices[0] + 1) * (num_slices[1] + 1);
          }
        }
      }; //<3, 2, 0>

      /// vert-at-hexadron for hexahedron mesh
      template<>
      struct StructIndexMapping<3, 3, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          return (i / (num_slices[0] * num_slices[1]) + j / 4) * ((num_slices[0] + 1)*(num_slices[1] + 1))
                 + ((i % (num_slices[0] * num_slices[1])) / num_slices[0] + (j % 4) / 2) * (num_slices[0] + 1)
                 + (i % (num_slices[0] * num_slices[1])) % num_slices[0] + j % 2;
        }
      }; //<3, 3, 0>

      /// edge-at-quad for hexahedron mesh
      template<>
      struct StructIndexMapping<3, 2, 1>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // auxiliary indices
          Index pos, node_0, node_0_x, node_0_y, node_0_z;

          // xy plain
          if( i < num_slices[0] * num_slices[1] * (num_slices[2] + 1))
          {
            // node 0
            node_0 = i + (i / num_slices[0]) + (i / (num_slices[0] * num_slices[1])) * (num_slices[0] + 1);

            if (j/2 == 0) //x direction edge
            {
              return node_0 - (node_0 / (num_slices[0] + 1)) + (j % 2) * num_slices[0];
            }
            else // y direction edge
            {
              // mesh coordinates
              node_0_x = node_0 % (num_slices[0] + 1);
              node_0_y = (node_0 / (num_slices[0] + 1)) % (num_slices[1] + 1);
              node_0_z = node_0 / ((num_slices[0] + 1) * (num_slices[1] + 1));
              return num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) + node_0_y + (node_0_x + j - 2)
                     * num_slices[1] + node_0_z * (num_slices[0] + 1) * num_slices[1];
            }
          }
          // xz plain
          else if( i >= num_slices[0] * num_slices[1] * (num_slices[2]+1)
                   && i < num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1)))
          {
            // position and node 0
            pos = i - num_slices[0] * num_slices[1] * (num_slices[2] + 1);
            node_0 = pos % num_slices[0] + (pos / (num_slices[0] * num_slices[2])) * (num_slices[0] + 1)
                     + (pos % (num_slices[0] * num_slices[2])) / num_slices[0]
                     * (num_slices[0] + 1) * (num_slices[1] + 1);

            if (j/2 == 0) //x direction edge
            {
              return node_0 - (node_0 / (num_slices[0] + 1)) + (j % 2) * num_slices[0] * (num_slices[1] + 1);
            }
            else // z direction edge
            {
              // coordinates
              node_0_x = node_0 % (num_slices[0] + 1);
              node_0_y = (node_0 / (num_slices[0] + 1)) % (num_slices[1] + 1);
              node_0_z = node_0 / ((num_slices[0] + 1) * (num_slices[1] + 1));
              return num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1)
                     + num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1)
                     + node_0_z + (node_0_x + j - 2) * num_slices[2] + node_0_y * (num_slices[0] + 1) * num_slices[2];
            }
          }
          // yz plain
          else
          {
            // position and node 0
            pos = i - num_slices[0] * (num_slices[1] * (num_slices[2] + 1) + num_slices[2] * (num_slices[1] + 1));
            node_0 = pos / (num_slices[1] * num_slices[2]) + (pos % num_slices[1]) * (num_slices[0] + 1)
                     + ((pos / num_slices[1]) % num_slices[2]) * (num_slices[0] + 1) * (num_slices[1] + 1);

            // coordinates
            node_0_x = node_0 % (num_slices[0] + 1);
            node_0_y = (node_0 / (num_slices[0] + 1)) % (num_slices[1] + 1);
            node_0_z = node_0 / ((num_slices[0] + 1) * (num_slices[1] + 1));

            if (j/2 == 0) //y direction edge
            {
              return num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) + node_0_y + node_0_x * num_slices[1]
                     + (node_0_z + j) * (num_slices[0] + 1) * num_slices[1];
            }
            else // z direction edge
            {
              return num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) + num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1)
                     + node_0_z + node_0_x * num_slices[2] + (node_0_y + j - 2) * (num_slices[0] + 1) * num_slices[2];
            }
          }//if

        }
      }; //<3, 2, 1>


      /// edge-at-hexadron for hexahedron mesh
      template<>
      struct StructIndexMapping<3, 3, 1>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          Index x = i % (num_slices[0] * num_slices[1]) % num_slices[0];
          Index y = i % (num_slices[0] * num_slices[1]) / num_slices[0];
          Index z = i / (num_slices[0] * num_slices[1]);
          if( j / 4 == 0)
          {
            return (z + j / 2) * num_slices[0] * (num_slices[1] + 1) + (y + j % 2) * num_slices[0] + x;
          }
          else if( j / 4 == 1)
          {
            return (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0] + (z + j / 6) * num_slices[1] * (num_slices[0]+1)
                    + (x + j % 2) * num_slices[1] + y;
          }
          else
          {
            return (num_slices[2] + 1) * (num_slices[1] + 1) * num_slices[0] + num_slices[1] * (num_slices[0] + 1)
                   * (num_slices[2] + 1) + (y + j / 10) * num_slices[2] * (num_slices[0] + 1) + (x + j % 2) * num_slices[2] + z;
          }
        }
      }; //<3, 3, 1>

      /// quad-at-hexadron for hexahedron mesh
      template<>
      struct StructIndexMapping<3, 3, 2>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {

          if(j/2 == 0) // x surface
          {
            return i + j*num_slices[0] * num_slices[1];
          }
          else if( j/4 == 0) // y surface
          {
            // compute index
            Index x = i % (num_slices[0] * num_slices[1]) % num_slices[0];
            Index y = i % (num_slices[0] * num_slices[1]) / num_slices[0];
            Index z = i / (num_slices[0] * num_slices[1]);
            return num_slices[0] * num_slices[1] * (num_slices[2] + 1)
                   + x + z * num_slices[0] + (y + j - 2) * num_slices[0] * num_slices[2];
          }
          else // z surface
          {
            // compute index
            Index x = i % (num_slices[0] * num_slices[1]) % num_slices[0];
            Index y = i % (num_slices[0] * num_slices[1]) / num_slices[0];
            Index z = i / (num_slices[0] * num_slices[1]);
            return num_slices[0] * (num_slices[1] * (num_slices[2]+1) + num_slices[2] * (num_slices[1]+1))
                   + y + z * num_slices[1] + (x + j - 4) * num_slices[1] * num_slices[2];
          }
        }
      }; //<3, 3, 2>

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP
