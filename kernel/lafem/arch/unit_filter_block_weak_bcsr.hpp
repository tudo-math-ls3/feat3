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
  struct UnitFilterBlockWeakBCSR
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool ign_nans, const int block_height, const int block_width);

    template<typename DT_, typename IT_>
    static void exec_cuda_impl(DT_* mat_a, const DT_* const mat_m, const IT_* const row_ptr, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool ign_nans, const int block_height, const int block_width);

    template<typename DT_, typename IT_, int block_height_, int block_width_>
    static void exec_generic(Memory::Arbiter& mat_a, const Memory::Arbiter& mat_m, const Memory::Arbiter& row_ptr, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool ign_nans)
    {
      Memory::TypedView<DT_> a_view(mat_a.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<DT_> m_view(mat_m.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> f_view(f_val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(a_view.get_w(), m_view.get_r(), rp_view.get_r(), f_view.get_r(), fi_view.get_r(), f_nzes, ign_nans, block_height_, block_width_);
    }

    template<typename DT_, typename IT_, int block_height_, int block_width_>
    static void exec_cuda(Memory::Arbiter& mat_a, const Memory::Arbiter& mat_m, const Memory::Arbiter& row_ptr, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool ign_nans)
    {
      Memory::TypedView<DT_> a_view(mat_a.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<DT_> m_view(mat_m.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> f_view(f_val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(a_view.get_w(), m_view.get_r(), rp_view.get_r(), f_view.get_r(), fi_view.get_r(), f_nzes, ign_nans, block_height_, block_width_);
    }

    template<typename DT_, typename IT_, int block_height_, int block_width_>
    static void exec(Memory::Arbiter& mat_a, const Memory::Arbiter& mat_m, const Memory::Arbiter& row_ptr, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool ign_nans)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_height_, block_width_>(mat_a, mat_m, row_ptr, f_val, f_idx, f_nzes, ign_nans);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::UnitFilterBlockWeakBCSR::exec");
        exec_generic<DT_, IT_, block_height_, block_width_>(mat_a, mat_m, row_ptr, f_val, f_idx, f_nzes, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockWeakBCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, block_height_, block_width_>(mat_a, mat_m, row_ptr, f_val, f_idx, f_nzes, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockWeakBCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::UnitFilterBlockWeakBCSR::exec: unknown backend!");
        break;
      }
    }
  }; // struct UnitFilterBlockWeakBCSR

#ifdef FEAT_EICKT
  extern template void UnitFilterBlockWeakBCSR::exec_generic_impl(float* mat_a, const float* mat_m, const std::uint32_t* const row_ptr, const float * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  extern template void UnitFilterBlockWeakBCSR::exec_generic_impl(float* mat_a,const float* mat_m, const std::uint64_t* const row_ptr, const float * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  extern template void UnitFilterBlockWeakBCSR::exec_generic_impl(double* mat_a, const double* mat_m, const std::uint32_t* const row_ptr, const double * const f_val, const std::uint32_t * const f_idx, const Index f_nzes, const bool, const int, const int);
  extern template void UnitFilterBlockWeakBCSR::exec_generic_impl(double* mat_a,const double* mat_m, const std::uint64_t* const row_ptr, const double * const f_val, const std::uint64_t * const f_idx, const Index f_nzes, const bool, const int, const int);
#endif
} // namespace FEAT::LAFEM::Arch
  /// \endcond

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCK_WEAK_BCSR_HPP 1
#include <kernel/lafem/arch/unit_filter_block_weak_bcsr_generic.hpp>
#endif
