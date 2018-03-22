// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template class ConformalMesh<Shape::Simplex<2>, 2, Real>;
    template class ConformalMesh<Shape::Simplex<3>, 3, Real>;
    template class ConformalMesh<Shape::Hypercube<2>, 2, Real>;
    template class ConformalMesh<Shape::Hypercube<3>, 3, Real>;

    template class StandardRefinery<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    template class StandardRefinery<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    template class StandardRefinery<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    template class StandardRefinery<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;
  } // namespace Geometry
} // namespace FEAT
