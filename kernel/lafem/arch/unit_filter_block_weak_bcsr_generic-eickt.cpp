// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/unit_filter_block_weak_bcsr.hpp>

namespace FEAT::LAFEM::Arch
{
  template void UnitFilterBlockWeakBCSR::exec_generic_impl(float* mat_a, const float* mat_m, const std::uint32_t* const row_ptr, const float * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  template void UnitFilterBlockWeakBCSR::exec_generic_impl(float* mat_a,const float* mat_m, const std::uint64_t* const row_ptr, const float * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  template void UnitFilterBlockWeakBCSR::exec_generic_impl(double* mat_a, const double* mat_m, const std::uint32_t* const row_ptr, const double * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  template void UnitFilterBlockWeakBCSR::exec_generic_impl(double* mat_a,const double* mat_m, const std::uint64_t* const row_ptr, const double * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool, const int, const int);
} // namespace FEAT::LAFEM::Arch
