#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_functional.hpp>

#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>

//#include <kernel/meshopt/rumpf_functionals/2d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_q1_unrolled.hpp>
//
//#include <kernel/meshopt/rumpf_functionals/3d_p1_unrolled.hpp>
//#include <kernel/meshopt/rumpf_functionals/3d_p1_d2.hpp>
//
//#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
//#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>

using namespace FEAT;

// RumpfFunctionals
template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>, 3, 3, double>>
>;

//template class FEAT::Meshopt::RumpfFunctionalUnrolled
//<
//  double,
//  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>
//>;
//
//template class FEAT::Meshopt::RumpfFunctionalUnrolled
//<
//  double,
//  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, 3, double>>
//>;
//
//template class FEAT::Meshopt::RumpfFunctionalUnrolled
//<
//  double,
//  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>
//>;

// Split variants
/// \compilerhack icc < 16.0 seems to be unable to deal with template template parameters and extern
//#if !(defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600))
//template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional>;
//template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional_D2>;
//
//template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional>;
//template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional_D2>;
//#endif
