// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>

#include <kernel/meshopt/rumpf_trafo.hpp>

using namespace FEAT;

// RumpfTrafos
template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<3>, 3, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, double>>, double>;

template struct FEAT::Meshopt::RumpfTrafo<
  Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<3>, 3, double>>, double>;
