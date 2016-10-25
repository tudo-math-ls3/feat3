#include <kernel/geometry/mesh_part.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template class MeshPart<ConformalMesh<Shape::Simplex<2>, 2, 2, Real>>;
    template class MeshPart<ConformalMesh<Shape::Simplex<3>, 3, 3, Real>>;
    template class MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, 2, Real>>;
    template class MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, 3, Real>>;

    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<2>, 2, 2, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<3>, 3, 3, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, 2, Real>>>;
    template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, 3, Real>>>;
  } // namespace Geometry
} // namespace FEAT
