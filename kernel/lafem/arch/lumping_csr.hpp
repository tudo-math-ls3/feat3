// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_arbiter.hpp>

namespace FEAT::LAFEM::Arch
{
  /// for 0 <= i < n: lump[i] <- sum_{0 <= j < n} A[i][j]
  struct LumpingCSR
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * lump, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows);

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * lump, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& lump, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows)
    {
      Memory::TypedView<DT_> l_view(lump.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(l_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows);
    }

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& lump, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows)
    {
      Memory::TypedView<DT_> l_view(lump.view(Memory::Location::cuda, Memory::Access::write));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(l_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows);
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& lump, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(lump, val, col_idx, row_ptr, rows);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::LumpingCSR::exec");
        exec_generic<DT_, IT_>(lump, val, col_idx, row_ptr, rows);
#else
        XABORTM("LAFEM::Arch::LumpingCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(lump, val, col_idx, row_ptr, rows);
#else
        XABORTM("LAFEM::Arch::LumpingCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::LumpingCSR::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template void LumpingCSR::exec_generic_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
  extern template void LumpingCSR::exec_generic_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
  extern template void LumpingCSR::exec_generic_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index);
  extern template void LumpingCSR::exec_generic_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_LUMPING_CSR_HPP 1
#include <kernel/lafem/arch/lumping_csr_generic.hpp>
#endif
