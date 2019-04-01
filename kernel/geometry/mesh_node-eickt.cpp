// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/geometry/mesh_node.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template class MeshNode<ConformalMesh<Shape::Simplex<2>, 2, Real>, ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    template class MeshNode<ConformalMesh<Shape::Simplex<3>, 3, Real>, ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    template class MeshNode<ConformalMesh<Shape::Hypercube<2>, 2, Real>, ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    template class MeshNode<ConformalMesh<Shape::Hypercube<3>, 3, Real>, ConformalMesh<Shape::Hypercube<3>, 3, Real>>;

    template class MeshNode<ConformalMesh<Shape::Simplex<2>, 2, Real>, MeshPart<ConformalMesh<Shape::Simplex<2>, 2, Real>>>;
    template class MeshNode<ConformalMesh<Shape::Simplex<3>, 3, Real>, MeshPart<ConformalMesh<Shape::Simplex<3>, 3, Real>>>;
    template class MeshNode<ConformalMesh<Shape::Hypercube<2>, 2, Real>, MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, Real>>>;
    template class MeshNode<ConformalMesh<Shape::Hypercube<3>, 3, Real>, MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, Real>>>;

    template class RootMeshNode<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    template class RootMeshNode<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    template class RootMeshNode<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    template class RootMeshNode<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;
  } // namespace Geometry
} // namespace FEAT
