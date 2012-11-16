#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP 1

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
      typedef ConformalMesh<Shape::Hexahedron> HexaMesh;

      HexaMesh* create_hexa_mesh_3d(int orientation);

      void validate_refined_hexa_mesh_3d(const HexaMesh& mesh, int orientation);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP
