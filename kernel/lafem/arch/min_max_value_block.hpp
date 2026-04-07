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
  struct MinMaxValueBlock
  {
    template <typename ValueType_>
    static ValueType_ exec_generic_impl(const ValueType_ * const x, const Index size, const bool min_, const bool abs_);

    template <typename ValueType_>
    static ValueType_ exec_generic(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      Memory::TypedView<ValueType_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(x_view.get_r(), size, min_, abs_);
    }

    template <typename ValueType_>
    static ValueType_ exec(const Memory::Arbiter& x, const Index size, const bool min_, const bool abs_)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<ValueType_>(x, size, min_, abs_);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MinMaxValueBlock::exec");
        return exec_generic<ValueType_>(x, size, min_, abs_);
#else
        XABORTM("LAFEM::Arch::MinMaxValueBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::MinMaxValueBlock::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::MinMaxValueBlock::exec: unknown backend!");
        break;
      }
    }
  }; // struct MinMaxValueBlock
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MIN_MAX_VALUE_BLOCK_HPP 1
#include <kernel/lafem/arch/min_max_value_block_generic.hpp>
#endif
