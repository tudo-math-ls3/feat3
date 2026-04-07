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
  /// for 0 <= i < n, 0 <= k < blocksize: r[i][k] <- r[i][k] + alpha[k]*x[i][k]
  struct AxpyBlock
  {
    template <typename ValueType_>
    static void exec_generic_impl(ValueType_ * r, const ValueType_ alpha, const ValueType_ * const x, const Index size);

    template<typename ValueType_>
    static void exec_generic(Memory::Arbiter& r, const ValueType_ alpha, const Memory::Arbiter& x, const Index size)
    {
      Memory::TypedView<ValueType_> r_view(r.view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap));
      Memory::TypedView<ValueType_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), size);
    }

    template<typename ValueType_>
    static void exec(Memory::Arbiter& r, const ValueType_ alpha, const Memory::Arbiter& x, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<ValueType_>(r, alpha, x, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::AxpyBlock::exec");
        exec_generic<ValueType_>(r, alpha, x, size);
#else
        XABORTM("LAFEM::Arch::AxpyBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::AxpyBlock::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::AxpyBlock::exec: unknown backend!");
        break;
      }
    }
  }; // struct AxpyBlock
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_AXPY_BLOCK_HPP 1
#include <kernel/lafem/arch/axpy_block_generic.hpp>
#endif
