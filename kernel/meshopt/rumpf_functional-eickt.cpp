// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_functional.hpp>

#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>

using namespace FEAT;

// RumpfFunctionals
template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<3>, 3, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, double>>
>;

template class FEAT::Meshopt::RumpfFunctional
<
  double,
  Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>, 3, double>>
>;
