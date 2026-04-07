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
  /// for all 0 <= k < blocksize: r[k] <- sum_{0 <= i < n} x[i][k] * y[i][k]
  struct DotProductBlock
  {
    template <typename ValueType_>
    static ValueType_ exec_generic_impl(const ValueType_ * const x, const ValueType_ * const y, const Index size);

    template <typename ValueType_>
    static ValueType_ exec_generic(const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      Memory::TypedView<ValueType_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<ValueType_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(x_view.get_r(), y_view.get_r(), size);
    }

    template <typename ValueType_>
    static ValueType_ exec(const Memory::Arbiter& x, const Memory::Arbiter& y, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<ValueType_>(x, y, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::DotProductBlock::exec");
        return exec_generic<ValueType_>(x, y, size);
#else
        XABORTM("LAFEM::Arch::DotProductBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::DotProductBlock::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::DotProductBlock::exec: unknown backend!");
        break;
      }
    }
  };
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_BLOCK_HPP 1
#include <kernel/lafem/arch/dot_product_block_generic.hpp>
#endif
