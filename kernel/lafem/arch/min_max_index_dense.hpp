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
  struct MinMaxIndexDense
  {
    template <typename DT_>
    static Index exec_generic_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_);

    static Index exec_mkl_impl(const float * const x, const Index size, const bool min_, const bool abs_);
    static Index exec_mkl_impl(const double * const x, const Index size, const bool min_, const bool abs_);

    template <typename DT_>
    static Index exec_mkl_impl(const DT_ * const, const Index, const bool, const bool)
    {
      XABORTM("LAFEM::Arch::MinMaxIndexDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static Index exec_cuda_impl(const DT_ * const x, const Index size, const bool min_, const bool abs_);

    template <typename DT_>
    static Index exec_generic(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(x_view.get_r(), size, min_, abs_);
    }

    template <typename DT_>
    static Index exec_mkl(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      return exec_mkl_impl(x_view.get_r(), size, min_, abs_);
    }

    template <typename DT_>
    static Index exec_cuda(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      return exec_cuda_impl(x_view.get_r(), size, min_, abs_);
    }

    template <typename DT_>
    static Index exec(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<DT_>(x, size, min_, abs_);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        return exec_mkl<DT_>(x, size, min_, abs_);
#else
        XABORTM("LAFEM::Arch::MinMaxIndexDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        return exec_cuda<DT_>(x, size, min_, abs_);
#else
        XABORTM("LAFEM::Arch::MinMaxIndexDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MinMaxIndexDense::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template Index MinMaxIndexDense::exec_generic_impl(const float * const, const Index, const bool, const bool);
  extern template Index MinMaxIndexDense::exec_generic_impl(const double * const, const Index, const bool, const bool);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MIN_MAX_INDEX_DENSE_HPP 1
#include <kernel/lafem/arch/min_max_index_dense_generic.hpp>
#endif
