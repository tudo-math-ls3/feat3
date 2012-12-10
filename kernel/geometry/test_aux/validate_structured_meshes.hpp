#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_VALIDATE_STRUCTURED_MESHES_HPP
#define KERNEL_GEOMETRY_TEST_AUX_VALIDATE_STRUCTURED_MESHES_HPP 1

// includes, FEAST
#include <kernel/geometry/structured_mesh.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {

      void validate_structured_mesh_2d(const StructuredMesh<2>* mesh);
      void validate_structured_mesh_3d(const StructuredMesh<3>* mesh);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_VALIDATE_STRUCTURED_MESHES_HPP