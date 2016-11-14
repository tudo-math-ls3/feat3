#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_q1_unrolled.hpp>

//#include <kernel/meshopt/rumpf_functionals/3d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>

#include <kernel/meshopt/hyperelasticity_functional.hpp>

using namespace FEAT;

// HyperelasticityFunctionals
template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
  Meshopt::RumpfFunctional
  <
    double,
    Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>
  >
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
  Meshopt::RumpfFunctional
  <
    double,
    Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>
  >
>;
