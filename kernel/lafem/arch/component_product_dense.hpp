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
  /// for 0 <= i < n: r[i] <- x[i] * y[i]
  struct ComponentProductDense
  {
    template <typename DT_>
    static void exec_generic_impl(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);

    static void exec_mkl_impl(float * r, const float * const x, const float * const y, const Index size);
    static void exec_mkl_impl(double * r, const double * const x, const double * const y, const Index size);

    template <typename DT_>
    static void exec_mkl_impl(DT_ *, const DT_ * const, const DT_ * const, const Index)
    {
      XABORTM("LAFEM::Arch::ComponentProductDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static void exec_cuda_impl(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);

    template <typename DT_>
    static void exec_generic(Memory::Arbiter& r, const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), x_view.get_r(), y_view.get_r(), size);
    }

    template <typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      exec_mkl_impl(r_view.get_w(), x_view.get_r(), y_view.get_r(), size);
    }

    template <typename DT_>
    static void exec_cuda(Memory::Arbiter& r, const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(r_view.get_w(), x_view.get_r(), y_view.get_r(), size);
    }

    template <typename DT_>
    static void exec(Memory::Arbiter& r, const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_>(r, x, y, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        exec_mkl<DT_>(r, x, y, size);
#else
        XABORTM("LAFEM::Arch::ComponentProduct::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_>(r, x, y, size);
#else
        XABORTM("LAFEM::Arch::ComponentProduct::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::ComponentProduct::exec: unknown backend!");
        break;
      }
    }
  }; // struct ComponentProductDense

#ifdef FEAT_EICKT
  extern template void ComponentProductDense::exec_generic_impl(float *, const float * const, const float * const, const Index);
  extern template void ComponentProductDense::exec_generic_impl(double *, const double * const, const double * const, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_COMPONENT_PRODUCT_DENSE_HPP 1
#include <kernel/lafem/arch/component_product_dense_generic.hpp>
#endif
