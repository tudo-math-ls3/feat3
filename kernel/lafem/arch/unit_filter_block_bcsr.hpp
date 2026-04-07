// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_arbiter.hpp>

/// \cond internal
namespace FEAT::LAFEM::Arch
{
  struct UnitFilterBlockBCSR
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);

    template<typename DT_, typename IT_, int bh_, int bw_>
    static void exec_generic(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const bool ign_nans)
    {
      Memory::TypedView<DT_> a_view(mat.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> fv_view(f_val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(a_view.get_w(), rp_view.get_r(), ci_view.get_r(), fv_view.get_r(), fi_view.get_r(), f_nzes, unit, ign_nans, bh_, bw_);
    }

    template<typename DT_, typename IT_>
    static void exec_cuda_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);

    template<typename DT_, typename IT_, int bh_, int bw_>
    static void exec_cuda(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const bool ign_nans)
    {
      Memory::TypedView<DT_> a_view(mat.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> fv_view(f_val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(a_view.get_w(), rp_view.get_r(), ci_view.get_r(), fv_view.get_r(), fi_view.get_r(), f_nzes, unit, ign_nans, bh_, bw_);
    }

    template<typename DT_, typename IT_, int bh_, int bw_>
    static void exec(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const bool ign_nans)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, bh_, bw_>(mat, row_ptr, col_idx, f_val, f_idx, f_nzes, unit, ign_nans);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::UnitFilterBlockBCSR::exec");
        exec_generic<DT_, IT_, bh_, bw_>(mat, row_ptr, col_idx, f_val, f_idx, f_nzes, unit, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockBCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, bh_, bw_>(mat, row_ptr, col_idx, f_val, f_idx, f_nzes, unit, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockBCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::UnitFilterBlockBCSR::exec: unknown backend!");
        break;
      }
    }
  }; // struct UnitFilterBlockBCSR

#ifdef FEAT_EICKT
  extern template void UnitFilterBlockBCSR::exec_generic_impl(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const float * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  extern template void UnitFilterBlockBCSR::exec_generic_impl(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const float * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  extern template void UnitFilterBlockBCSR::exec_generic_impl(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const double * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
  extern template void UnitFilterBlockBCSR::exec_generic_impl(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const double * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const bool ign_nans, const int block_height, const int block_width);
#endif
} // namespace FEAT::LAFEM::Arch
/// \endcond

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCK_BCSR_HPP 1
#include <kernel/lafem/arch/unit_filter_block_bcsr_generic.hpp>
#endif
