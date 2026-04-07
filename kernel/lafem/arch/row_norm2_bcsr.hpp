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
  struct RowNorm2BCSR
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_* row_norms, const DT_* const a, const IT_* const col_idx,
      const IT_ * const row_ptr, const Index rows, const int bh, const int bw, const bool squared);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& row_norms, const Memory::Arbiter& a, const Memory::Arbiter& col_idx,
      const Memory::Arbiter&row_ptr, const Index rows, const int bh, const int bw, const bool squared)
    {
      Memory::TypedView<DT_> r_view(row_norms.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows, bh, bw, squared);
    }

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_* row_norms, const DT_* const a, const IT_* const col_idx,
      const IT_* const row_ptr, const Index rows, const int bh, const int bw, const bool squared);

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& row_norms, const Memory::Arbiter& a, const Memory::Arbiter& col_idx,
      const Memory::Arbiter&row_ptr, const Index rows, const int bh, const int bw, const bool squared)
    {
      Memory::TypedView<DT_> r_view(row_norms.view(Memory::Location::cuda, Memory::Access::write));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(r_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows, bh, bw, squared);
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& row_norms, const Memory::Arbiter& a, const Memory::Arbiter& col_idx,
      const Memory::Arbiter& row_ptr, const Index rows, const int bh, const int bw, const bool squared)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(row_norms, a, col_idx, row_ptr, rows, bh, bw, squared);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::RowNorm2BCSR::exec");
        exec_generic<DT_, IT_>(row_norms, a, col_idx, row_ptr, rows, bh, bw, squared);
#else
        XABORTM("LAFEM::Arch::RowNorm2BCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(row_norms, a, col_idx, row_ptr, rows, bh, bw, squared);
#else
        XABORTM("LAFEM::Arch::RowNorm2BCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::RowNorm2BCSR::exec: unknown backend!");
        break;
      }
    }
  }; // struct RowNorm2BCSR

#ifdef FEAT_EICKT
  extern template void RowNorm2BCSR::exec_generic_impl(float*, const float* const, const Index* const, const Index * const, const Index, const int, const int, const bool);
  extern template void RowNorm2BCSR::exec_generic_impl(double*, const double* const, const Index* const, const Index* const, const Index, const int, const int, const bool);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_ROW_NORM2_BCSR_HPP 1
#include <kernel/lafem/arch/row_norm2_bcsr_generic.hpp>
#endif
