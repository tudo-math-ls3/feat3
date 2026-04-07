// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/dot_product_dense.hpp>

namespace FEAT::LAFEM::Arch
{
  template float DotProductDense::exec_generic_impl(const float * const, const float * const, const Index);
  template double DotProductDense::exec_generic_impl(const double * const, const double * const, const Index);
} // namespace FEAT::LAFEM::Arch
