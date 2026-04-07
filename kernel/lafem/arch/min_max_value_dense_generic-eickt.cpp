// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/min_max_value_dense.hpp>

namespace FEAT::LAFEM::Arch
{
  template float MinMaxValueDense::exec_generic_impl(const float * const, const Index, const bool, const bool);
  template double MinMaxValueDense::exec_generic_impl(const double * const, const Index, const bool, const bool);
} // namespace FEAT::LAFEM::Arch
