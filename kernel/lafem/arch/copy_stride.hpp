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
  /// for 0 <= i < n: r[i*stride_r + comp_r] <- x[i*stride_x + comp_x]
  struct CopyStride
  {
    template <typename DT_>
    static void exec_generic_impl(DT_ * r, const DT_ * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count);

    static void exec_mkl_impl(float * r, const float * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count);
    static void exec_mkl_impl(double * r, const double * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count);

    template <typename DT_>
    static void exec_mkl_impl(DT_ *, const DT_ * const, const int, const int, const int, const int, const Index)
    {
      XABORTM("LAFEM::Arch::CopyStride::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static void exec_cuda_impl(DT_ * r, const DT_ * const x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count);

    template <typename DT_>
    static void exec_generic(Memory::Arbiter& r, const Memory::Arbiter& x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, (stride_r > 1 ? Memory::Access::read_write : Memory::Access::write)));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), x_view.get_r(), stride_r, comp_r, stride_x, comp_x, count);
    }

    template <typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const Memory::Arbiter& x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, (stride_r > 1 ? Memory::Access::read_write : Memory::Access::write)));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      exec_mkl_impl(r_view.get_w(), x_view.get_r(), stride_r, comp_r, stride_x, comp_x, count);
    }

    template <typename DT_>
    static void exec_cuda(Memory::Arbiter& r, const Memory::Arbiter& x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, (stride_r > 1 ? Memory::Access::read_write : Memory::Access::write)));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(r_view.get_w(), x_view.get_r(), stride_r, comp_r, stride_x, comp_x, count);
    }

    template <typename DT_>
    static void exec(Memory::Arbiter& r, const Memory::Arbiter& x, const int stride_r, const int comp_r, const int stride_x, const int comp_x, const Index count)
    {
      XASSERT(r != x);
      XASSERT(comp_r < stride_r);
      XASSERT(comp_x < stride_x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_>(r, x, stride_r, comp_r, stride_x, comp_x, count);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        exec_mkl<DT_>(r, x, stride_r, comp_r, stride_x, comp_x, count);
#else
        XABORTM("LAFEM::Arch::CopyStride::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_>(r, x, stride_r, comp_r, stride_x, comp_x, count);
#else
        XABORTM("LAFEM::Arch::CopyStride::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::CopyStride::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template void CopyStride::exec_generic_impl(float *, const float * const, const int, const int, const int, const int, const Index);
  extern template void CopyStride::exec_generic_impl(double *, const double * const, const int, const int, const int, const int, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_COPY_STRIDE_HPP 1
#include <kernel/lafem/arch/copy_stride_generic.hpp>
#endif
