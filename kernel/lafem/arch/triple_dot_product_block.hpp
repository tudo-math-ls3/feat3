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
  struct TripleDotProductBlock
  {
    template <typename ValueType_>
    static ValueType_ exec_generic_impl(const ValueType_ * const x, const ValueType_ * const y, const ValueType_ * const z, const Index size);

    template <typename ValueType_>
    static ValueType_ exec_generic(const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index size)
    {
      Memory::TypedView<ValueType_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<ValueType_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<ValueType_> z_view(z.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(x_view.get_r(), y_view.get_r(), z_view.get_r(), size);
    }

    template <typename ValueType_>
    static ValueType_ exec(const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index size)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<ValueType_>(x, y, z, size);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::TripleDotProductBlock::exec");
        return exec_generic<ValueType_>(x, y, z, size);
#else
        XABORTM("LAFEM::Arch::TripleDotProductBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::TripleDotProductBlock::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::TripleDotProductBlock::exec: unknown backend!");
        break;
      }
    }
  }; // struct TripleDotProductBlock
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_TRIPLE_DOT_PRODUCT_BLOCK_HPP 1
#include <kernel/lafem/arch/triple_dot_product_block_generic.hpp>
#endif
