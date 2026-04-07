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
  /// for 0 <= i < n, 0 <= j < blocksize: (diag[i])[j] <- (A[i][i])[j][j]
  struct DiagonalBCSR
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_* diag, const DT_* const val, const IT_* const col_idx, const IT_ * const row_ptr, const Index rows, const int block_size);

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * diag, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const int block_soze);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& diag, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const int block_size)
    {
      Memory::TypedView<DT_> d_view(diag.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(d_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows, block_size);
    }

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& diag, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const int block_size)
    {
      Memory::TypedView<DT_> d_view(diag.view(Memory::Location::cuda, Memory::Access::write));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(d_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), rows, block_size);
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& diag, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const int block_size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(diag, val,  col_idx, row_ptr, rows, block_size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::DiagonalBCSR::exec");
        exec_generic<DT_, IT_>(diag, val,  col_idx, row_ptr, rows, block_size);
#else
        XABORTM("LAFEM::Arch::DiagonalBCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(diag, val,  col_idx, row_ptr, rows, block_size);
#else
        XABORTM("LAFEM::Arch::DiagonalBCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::DiagonalBCSR::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template void DiagonalBCSR::exec_generic_impl(float *, const float * const, const Index * const, const Index * const, const Index, const int);
  extern template void DiagonalBCSR::exec_generic_impl(double *, const double * const, const Index * const, const Index * const, const Index, const int);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_DIAGONAL_BCSR_HPP 1
#include <kernel/lafem/arch/diagonal_bcsr_generic.hpp>
#endif
