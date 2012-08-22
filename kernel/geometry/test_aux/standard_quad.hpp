#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      typedef ConformalMesh< ConformalMeshPolicy< Shape::Quadrilateral > > QuadMesh;
      typedef ConformalSubMesh< ConformalSubMeshPolicy< Shape::Quadrilateral > > QuadSubMesh;

      QuadMesh* create_quad_mesh_2d(int orientation);

      void validate_refined_quad_mesh_2d(const QuadMesh& mesh, int orientation);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP
