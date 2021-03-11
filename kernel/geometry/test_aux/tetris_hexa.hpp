// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      typedef ConformalMesh<Shape::Hexahedron> HexaMesh;
      typedef MeshPart<HexaMesh> HexaSubMesh;

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
          y                   |   /               |   /
          ^                   |  J                |  K
          |                   | /                 | /
          |                   |/                  |/
         11---->x             0---------A---------1
         /
       z/
       v
      * \endverbatim
      * \author Peter Zajac
      */
      HexaMesh* create_tetris_mesh_3d();

      void validate_refined_tetris_mesh_3d(const HexaMesh& mesh);

      HexaMesh* create_big_tetris_mesh_3d();

      HexaSubMesh* create_tetris_quad_submesh_3d();

      void validate_refined_tetris_quad_submesh_3d(const HexaSubMesh& mesh);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_TETRIS_HEXA_HPP
