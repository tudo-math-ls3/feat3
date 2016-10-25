#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAT
{
  namespace Cubature
  {
    template bool DynamicFactory::create<Shape::Simplex<2>, Real, Real, Tiny::Vector<Real, 2, 2>>
      (Rule<Shape::Simplex<2>, Real, Real, Tiny::Vector<Real, 2, 2>>&) const;
    template bool DynamicFactory::create<Shape::Simplex<3>, Real, Real, Tiny::Vector<Real, 3, 3>>
      (Rule<Shape::Simplex<3>, Real, Real, Tiny::Vector<Real, 3, 3>>&) const;
    template bool DynamicFactory::create<Shape::Hypercube<2>, Real, Real, Tiny::Vector<Real, 2, 2>>
      (Rule<Shape::Hypercube<2>, Real, Real, Tiny::Vector<Real, 2, 2>>&) const;
    template bool DynamicFactory::create<Shape::Hypercube<3>, Real, Real, Tiny::Vector<Real, 3, 3>>
      (Rule<Shape::Hypercube<3>, Real, Real, Tiny::Vector<Real, 3, 3>>&) const;
  } // namespace Cubature
} // namespace FEAT
