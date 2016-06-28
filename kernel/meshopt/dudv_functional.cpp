#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/dudv_functional.hpp>

using namespace FEAT;

// DuDvFunctionals
template class FEAT::Meshopt::DuDvFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>,
  LAFEM::SparseMatrixBCSR
>;

template class FEAT::Meshopt::DuDvFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
  LAFEM::SparseMatrixBCSR
>;
