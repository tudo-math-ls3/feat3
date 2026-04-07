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
  struct TripleDotProductDense
  {
    template <typename DT_>
    static DT_ exec_generic_impl(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);

    template <typename DT_>
    static DT_ exec_cuda_impl(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);

    template <typename DT_>
    static DT_ exec_generic(const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index size)
    {
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> z_view(z.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(x_view.get_r(), y_view.get_r(), z_view.get_r(), size);
    }

    template <typename DT_>
    static DT_ exec_cuda(const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index size)
    {
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> z_view(z.view(Memory::Location::cuda, Memory::Access::read));
      return exec_cuda_impl(x_view.get_r(), y_view.get_r(), z_view.get_r(), size);
    }

    template <typename DT_>
    static DT_ exec(const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<DT_>(x, y, z, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::TripleDotProductDense::exec");
        return exec_generic<DT_>(x, y, z, size);
#else
        XABORTM("LAFEM::Arch::TripleDotProductDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        return exec_cuda<DT_>(x, y, z, size);
#else
        XABORTM("LAFEM::Arch::TripleDotProductDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::TripleDotProductDense::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template float TripleDotProductDense::exec_generic_impl(const float * const, const float * const, const float * const, const Index);
  extern template double TripleDotProductDense::exec_generic_impl(const double * const, const double * const, const double * const, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_TRIPLE_DOT_PRODUCT_DENSE_HPP 1
#include <kernel/lafem/arch/triple_dot_product_dense_generic.hpp>
#endif
