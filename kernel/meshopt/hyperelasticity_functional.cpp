#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>

#include <kernel/meshopt/hyperelasticity_functional.hpp>

using namespace FEAT;

// HyperelasticityFunctionals
template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
  Meshopt::RumpfFunctional<double,Shape::Simplex<2>>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
  Meshopt::RumpfFunctional<double,Shape::Hypercube<2>>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>,
  Meshopt::RumpfFunctional_D2<double,Shape::Simplex<2>>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
  Meshopt::RumpfFunctional_D2<double,Shape::Hypercube<2>>
>;

// HyperelasticityFunctionals, local split functionals
template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
  Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
  Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional_D2>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>,
  Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional>
>;

template class FEAT::Meshopt::HyperelasticityFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>,
  Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional_D2>
>;
