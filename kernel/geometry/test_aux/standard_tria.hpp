// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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

      typedef ConformalMesh<Shape::Triangle> TriaMesh;
      typedef MeshPart<TriaMesh> TriaSubMesh;

      TriaMesh* create_tria_mesh_2d(int orientation);

      void validate_refined_tria_mesh_2d(const TriaMesh& mesh, int orientation);

      TriaMesh* create_patch_tria_mesh_2d();

      void validate_refined_patch_tria_mesh_2d(const TriaMesh& mesh);

      TriaSubMesh* create_patch_tria_submesh_2d();

      void validate_refined_patch_tria_submesh_2d(const TriaSubMesh& mesh);

      TriaSubMesh* create_patch_edge_submesh_2d();

      void validate_refined_patch_edge_submesh_2d(const TriaSubMesh& mesh);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
