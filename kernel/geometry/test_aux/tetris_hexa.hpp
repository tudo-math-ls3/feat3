#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      typedef ConformalMesh< ConformalMeshPolicy< Shape::Hexahedron > > HexaMesh;
      typedef ConformalSubMesh< ConformalSubMeshPolicy< Shape::Hexahedron > > HexaSubMesh;

      /**
       * \brief Creates the 3D tetris mesh.
       *
       * \verbatim
               13---------Z--------14---------A'-------15
               /|                  /|                  /|
              / |                 / |    18---------I'/-|----19
             O  |                P  |    /|          Q  |    /|
            /   W               /   X   / |         /   Y   / |
           /    |              /    |  D' |        /    |  E' |
          5---------H---------6---------I---------7     | /   H'
          |     |             |     |/    G'      |     |/    |
          |    10---------U---|----11---------V---|----12     |
          |    /              |    /|     |       |    /|     |
          E   /               F   / |    16-------G-F'/-|----17
          |  L                |  M  |    /        |  N  |    /
          | /                 | /   S   /         | /   T   /
          |/                  |/    |  B'         |/    |  C'
          2---------C---------3---------D---------4     | /
                              |     |/            |     |/
                              |     8---------R---|-----9
                              |    /              |    /
                              |   /               |   /
                              |  J                |  K
                              | /                 | /
                              |/                  |/
                              0---------A---------1 \endverbatim
      *
      * \author Peter Zajac
      */
      HexaMesh* create_tetris_mesh_3d()
      {
        Index num_entities[] =
        {
          20, // vertices
          36, // edges
          21, // quads
           4  // hexas
        };

        // create mesh
        HexaMesh* mesh = new HexaMesh(num_entities);

        // set up vertex coordinates array
        Real vtx[20*3] =
        {
           0.0, -1.0,  1.0,
           1.0, -1.0,  1.0,
          -1.0,  0.0,  1.0,
           0.0,  0.0,  1.0,
           1.0,  0.0,  1.0,
          -1.0,  1.0,  1.0,
           0.0,  1.0,  1.0,
           1.0,  1.0,  1.0,
           0.0, -1.0,  0.0,
           1.0, -1.0,  0.0,
          -1.0,  0.0,  0.0,
           0.0,  0.0,  0.0,
           1.0,  0.0,  0.0,
          -1.0,  1.0,  0.0,
           0.0,  1.0,  0.0,
           1.0,  1.0,  0.0,
           0.0, -1.0, -1.0,
           1.0, -1.0, -1.0,
           0.0,  0.0, -1.0,
           1.0,  0.0, -1.0
        };
        copy_vtx(mesh->get_vertex_set(), vtx);

        // set up vertices-at-edge
        Index v_e[36*2] =
        {
           0,  1,
           0,  3,
           1,  4,
           2,  3,
           3,  4,
           2,  5,
           3,  6,
           4,  7,
           5,  6,
           6,  7,
           8,  0,
           9,  1,
          10,  2,
          11,  3,
          12,  4,
          13,  5,
          14,  6,
          15,  7,
           8,  9,
           8, 11,
           9, 12,
          10, 11,
          11, 12,
          10, 13,
          11, 14,
          12, 15,
          13, 14,
          14, 15,
          16,  8,
          17,  9,
          18, 11,
          19, 12,
          16, 17,
          16, 18,
          17, 19,
          18, 19,
        };
        copy_idx(mesh->get_index_set<1,0>(), v_e);

        // set up vertices-at-face
        Index v_f[21*4] =
        {
           0,  1,  3,  4,
           2,  3,  5,  6,
           3,  4,  6,  7,
           8,  9,  0,  1,
           8, 11,  0,  3,
           9, 14,  1,  4,
          10, 11,  2,  3,
          11, 12,  3,  4,
          10, 13,  2,  5,
          11, 14,  3,  6,
          12, 15,  4,  7,
          13, 14,  5,  6,
          14, 15,  6,  7,
           8 , 9, 11, 12,
          10, 11, 13, 14,
          11, 12, 14, 15,
          16, 17,  8,  9,
          16, 18,  8, 11,
          17, 19,  9, 12,
          18, 19, 11, 12,
          16, 17, 18, 19
        };
        copy_idx(mesh->get_index_set<2,0>(), v_f);

        // set up edges-at-face
        Index e_f[21*4] =
        {
           0,  4,  1,  2,
           3,  8,  5,  6,
           4,  9,  6,  7,
          18,  0, 10, 11,
          19,  1, 10, 13,
          20,  2, 11, 14,
          21,  3, 12, 13,
          22,  4, 13, 14,
          23,  5, 12, 15,
          24,  6, 13, 17,
          25,  7, 14, 17,
          26,  8, 15, 16,
          27,  9, 16, 17,
          18, 22, 19, 20,
          21, 26, 23, 24,
          22, 27, 24, 25,
          32, 18, 28, 29,
          35, 22, 30, 31,
          33, 19, 28, 30,
          34, 20, 29, 31,
          32, 35, 33, 34
        };
        copy_idx(mesh->get_index_set<2,1>(), e_f);

        // set up vertices-at-cell
        Index v_c[4*8] =
        {
           8,  9, 11, 12,  0,  1,  3,  4,
          10, 11, 13, 14,  2,  3,  5,  6,
          11, 12, 14, 15,  3,  4,  6,  7,
          16, 17, 18, 19,  8,  9, 11, 12
        };
        copy_idx(mesh->get_index_set<3, 0>(), v_c);

        // set up edges-at-cell
        Index e_c[4*12] =
        {
          18, 22,  0,  4, 19, 20,  1,  2, 10, 11, 13, 14,
          21, 26,  3,  8, 23, 24,  5,  6, 12, 13, 15, 16,
          22, 27,  4,  9, 24, 25,  6,  7, 13, 14, 16, 17,
          32, 35, 18, 22, 33, 34, 19, 20, 28, 29, 30, 31
        };
        copy_idx(mesh->get_index_set<3,1>(), e_c);

        // set up faces-at-cell
        Index f_c[4*6] =
        {
          13,  0,  3,  7,  4,  5,
          14,  1,  6, 11,  8,  9,
          15,  2,  7, 12,  9, 10,
          20, 13, 16, 19, 17, 18
        };
        copy_idx(mesh->get_index_set<3,2>(), f_c);

        // okay
        return mesh;
      }
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP
