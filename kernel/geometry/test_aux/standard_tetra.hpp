#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {

      typedef ConformalMesh< ConformalMeshPolicy< Shape::Tetrahedron > > TetraMesh;
      typedef ConformalSubMesh< ConformalSubMeshPolicy< Shape::Tetrahedron > > TetraSubMesh;

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
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_TETRA_HPP