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
  struct ScaleColBCSR
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * r, const DT_ * const a, const IT_ * const col_idx, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index cols, const Index nonzeros, const int bh, const int bw);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& r, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Memory::Arbiter& x, const Index rows, const Index cols, const Index nonzeros, const int bh, const int bw)
    {
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), x_view.get_r(), rows, cols, nonzeros, bh, bw);
    }

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * r, const DT_ * const a, const IT_ * const col_idx, const IT_ * const row_ptr, const DT_ * const x, const Index rows, const Index, const Index, const int bh, const int bw);

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& r, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Memory::Arbiter& x, const Index rows, const Index cols, const Index nonzeros, const int bh, const int bw)
    {
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(r_view.get_w(), a_view.get_r(), ci_view.get_r(), rp_view.get_r(), x_view.get_r(), rows, cols, nonzeros, bh, bw);
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& r, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Memory::Arbiter& x, const Index rows, const Index cols, const Index nonzeros, const int bh, const int bw)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(r, a, col_idx, row_ptr, x, rows, cols, nonzeros, bh, bw);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::ScaleColBCSR::exec");
        exec_generic<DT_, IT_>(r, a, col_idx, row_ptr, x, rows, cols, nonzeros, bh, bw);
#else
        XABORTM("LAFEM::Arch::ScaleColBCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(r, a, col_idx, row_ptr, x, rows, cols, nonzeros, bh, bw);
#else
        XABORTM("LAFEM::Arch::ScaleColBCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::ScaleColBCSR::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template void ScaleColBCSR::exec_generic_impl(float *, const float * const, const std::uint64_t * const, const std::uint64_t * const, const float * const, const Index, const Index, const Index, const int, const int);
  extern template void ScaleColBCSR::exec_generic_impl(double *, const double * const, const std::uint64_t * const, const std::uint64_t * const, const double * const, const Index, const Index, const Index, const int, const int);
  extern template void ScaleColBCSR::exec_generic_impl(float *, const float * const, const std::uint32_t * const, const std::uint32_t * const, const float * const, const Index, const Index, const Index, const int, const int);
  extern template void ScaleColBCSR::exec_generic_impl(double *, const double * const, const std::uint32_t * const, const std::uint32_t * const, const double * const, const Index, const Index, const Index, const int, const int);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_SCALE_COL_BCSR_HPP 1
#include <kernel/lafem/arch/scale_col_bcsr_generic.hpp>
#endif
