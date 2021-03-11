// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAT
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
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_HEXA_HPP
