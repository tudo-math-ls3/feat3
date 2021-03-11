// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP 1

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

      typedef ConformalMesh<Shape::Tetrahedron> TetraMesh;
      typedef MeshPart<TetraMesh> TetraSubMesh;

      TetraMesh* create_tetra_mesh_3d(int orientation);

      void validate_refined_tetra_mesh_3d(const TetraMesh& mesh, int orientation);

      TetraMesh* create_big_tetra_mesh_3d();

      void validate_refined_big_tetra_mesh_3d(const TetraMesh& mesh);

      TetraMesh* create_really_big_tetra_mesh_3d();

      TetraSubMesh* create_tria_submesh_3d();

      void validate_refined_tria_submesh_3d(const TetraSubMesh& mesh);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP
