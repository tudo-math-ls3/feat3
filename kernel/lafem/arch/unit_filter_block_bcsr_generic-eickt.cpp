// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_block_bcsr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void UnitFilterBlockBCSR::exec_generic_impl(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const float * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_generic_impl(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const float * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_generic_impl(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const double * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  template void UnitFilterBlockBCSR::exec_generic_impl(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const double * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
} // namespace FEAT::LAFEM::Arch
