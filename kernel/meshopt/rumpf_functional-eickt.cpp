#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>

using namespace FEAT;

// RumpfFunctionals
template class FEAT::Meshopt::RumpfFunctional<double,Shape::Simplex<2>>;
template class FEAT::Meshopt::RumpfFunctional<double,Shape::Hypercube<2>>;

template class FEAT::Meshopt::RumpfFunctional_D2<double,Shape::Simplex<2>>;
template class FEAT::Meshopt::RumpfFunctional_D2<double,Shape::Hypercube<2>>;

// Split variants
/// \compilerhack icc < 16.0 seems to be unable to deal with template template parameters and extern
#if !(defined(FEAT_COMPILER_INTEL) && (FEAT_COMPILER_INTEL < 1600))
template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional>;
template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional_D2>;

template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional>;
template class FEAT::Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional_D2>;
#endif
