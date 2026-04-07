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
  /// for 0 <= i < n: r[i] <- r[i] + alpha*x[i]
  struct AxpyDense
  {
    template <typename DT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const Index size);

    static void exec_mkl_impl(float * r, const float alpha, const float * const x, const Index size);
    static void exec_mkl_impl(double * r, const double alpha, const double * const x, const Index size);

    template <typename DT_>
    static void exec_mkl_impl(DT_ *, const DT_, const DT_ * const, const Index)
    {
      XABORTM("LAFEM::Arch::AxpyDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const Index size);

    template<typename DT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), size);
    }

    template<typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), size);
    }

    template<typename DT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::read_write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), size);
    }

    template<typename DT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_>(r, alpha, x, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        exec_mkl<DT_>(r, alpha, x, size);
#else
        XABORTM("LAFEM::Arch::AxpyDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_>(r, alpha, x, size);
#else
        XABORTM("LAFEM::Arch::AxpyDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::AxpyDense::exec: unknown backend!");
        break;
      }
    }
  }; // struct AxpyDense

#ifdef FEAT_EICKT
  extern template void AxpyDense::exec_generic_impl(float *, const float, const float * const, const Index);
  extern template void AxpyDense::exec_generic_impl(double *, const double, const double * const, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_AXPY_DENSE_HPP 1
#include <kernel/lafem/arch/axpy_dense_generic.hpp>
#endif
