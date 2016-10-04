#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_trafo.hpp>

using namespace FEAT;

// RumpfTrafos
template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<3>, 3, 3, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<3>, 3, 3, double>>, double>;
