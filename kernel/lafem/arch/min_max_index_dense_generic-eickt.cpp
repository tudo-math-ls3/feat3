// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/min_max_index_dense.hpp>

namespace FEAT::LAFEM::Arch
{
  template Index MinMaxIndexDense::exec_generic_impl(const float * const, const Index, const bool, const bool);
  template Index MinMaxIndexDense::exec_generic_impl(const double * const, const Index, const bool, const bool);
} // namespace FEAT::LAFEM::Arch
