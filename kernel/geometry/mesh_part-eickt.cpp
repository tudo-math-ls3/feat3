// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/mesh_part.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template class MeshPart<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    template class MeshPart<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    template class MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    template class MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;

    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<2>, 2, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<3>, 3, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, Real>>>;
  } // namespace Geometry
} // namespace FEAT
