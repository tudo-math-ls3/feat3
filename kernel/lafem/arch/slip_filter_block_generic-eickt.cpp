// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/slip_filter_block.hpp>

namespace FEAT::LAFEM::Arch
{
  template void SlipFilterBlock::exec_generic_impl<float, std::uint64_t>(float * v, const float * const f_nu, const std::uint64_t * const f_idx, const Index f_nzes, const int bs);
  template void SlipFilterBlock::exec_generic_impl<double, std::uint64_t>(double * v, const double * const f_nu, const std::uint64_t * const f_idx, const Index f_nzes, const int bs);
  template void SlipFilterBlock::exec_generic_impl<float, std::uint32_t>(float * v, const float * const f_nu, const std::uint32_t * const f_idx, const Index f_nzes, const int bs);
  template void SlipFilterBlock::exec_generic_impl<double, std::uint32_t>(double * v, const double * const f_nu, const std::uint32_t * const f_idx, const Index f_nzes, const int bs);
} // namespace FEAT::LAFEM::Arch
