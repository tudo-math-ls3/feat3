#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_TRIA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_TRIA_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_sub_mesh.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {

      typedef ConformalMesh< ConformalMeshPolicy< Shape::Triangle > > TriaMesh;
      typedef ConformalSubMesh< ConformalSubMeshPolicy< Shape::Triangle > > TriaSubMesh;

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
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_TRIA_HPP
