// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/structured_mesh.hpp>

namespace FEAT
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
} // namespace FEAT
