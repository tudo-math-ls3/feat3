// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/dudv_functional.hpp>

using namespace FEAT;

// DuDvFunctionals
template class FEAT::Meshopt::DuDvFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, double>>,
  LAFEM::SparseMatrixBCSR
>;

template class FEAT::Meshopt::DuDvFunctional
<
  Mem::Main, double, Index,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, double>>,
  LAFEM::SparseMatrixBCSR
>;
