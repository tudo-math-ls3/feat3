// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/geometry/mesh_file_reader.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template class MeshNodeLinker<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    template class MeshNodeLinker<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    template class MeshNodeLinker<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    template class MeshNodeLinker<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;

    template void MeshFileReader::parse<ConformalMesh<Shape::Simplex<2>, 2, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Simplex<2>, 2, Real>>&,
      PartitionSet*);
    template void MeshFileReader::parse<ConformalMesh<Shape::Simplex<3>, 3, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Simplex<3>, 3, Real>>&,
      PartitionSet*);
    template void MeshFileReader::parse<ConformalMesh<Shape::Hypercube<2>, 2, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Hypercube<2>, 2, Real>>&,
      PartitionSet*);
    template void MeshFileReader::parse<ConformalMesh<Shape::Hypercube<3>, 3, Real>>(
      MeshNodeLinker<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      RootMeshNode<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      MeshAtlas<ConformalMesh<Shape::Hypercube<3>, 3, Real>>&,
      PartitionSet*);
  } // namespace Geometry
} // namespace FEAT
