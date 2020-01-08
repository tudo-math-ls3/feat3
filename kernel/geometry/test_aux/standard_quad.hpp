// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP
#define KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP 1

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
      typedef ConformalMesh<Shape::Quadrilateral> QuadMesh;
      typedef MeshPart<QuadMesh> QuadSubMesh;

      QuadMesh* create_quad_mesh_2d(int orientation);

      void validate_refined_quad_mesh_2d(const QuadMesh& mesh, int orientation);

    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_STANDARD_QUAD_HPP
