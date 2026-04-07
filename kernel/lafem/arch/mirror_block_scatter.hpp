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
  struct MirrorBlockScatter
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

    template<typename DT_, typename IT_>
    static void exec_cuda_impl(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

    template<typename DT_, typename IT_, int block_size_>
    static void exec_generic(const Index boff, const Index nidx, const Memory::Arbiter& idx, const Memory::Arbiter& buf, Memory::Arbiter& vec, const DT_ alpha)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(Index(block_size_), boff, nidx, i_view.get_r(), b_view.get_r(), x_view.get_w(), alpha);
    }

    template<typename DT_, typename IT_, int block_size_>
    static void exec_cuda(const Index boff, const Index nidx, const Memory::Arbiter& idx, const Memory::Arbiter& buf, Memory::Arbiter& vec, const DT_ alpha)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::cuda, Memory::Access::read));
      return exec_cuda_impl(Index(block_size_), boff, nidx, i_view.get_r(), b_view.get_r(), x_view.get_w(), alpha);
    }

    template<typename DT_, typename IT_, int block_size_>
    static void exec(const Index boff, const Index nidx, const Memory::Arbiter& idx, const Memory::Arbiter& buf, Memory::Arbiter& vec, const DT_ alpha)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_size_>(boff, nidx, idx, buf, vec, alpha);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MirrorBlockScatter::exec");
        exec_generic<DT_, IT_, block_size_>(boff, nidx, idx, buf, vec, alpha);
#else
        XABORTM("LAFEM::Arch::MirrorBlockScatter::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, block_size_>(boff, nidx, idx, buf, vec, alpha);
#else
        XABORTM("LAFEM::Arch::MirrorBlockScatter::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MirrorBlockScatter::exec: unknown backend!");
        break;
      }
    }
  }; // struct MirrorBlockScatter
} // namespace FEAT::LAFEM::Arch

#ifndef __CUDACC__
#define KERNEL_LAFEM_ARCH_MIRROR_BLOCK_SCATTER_HPP 1
#include <kernel/lafem/arch/mirror_block_scatter_generic.hpp>
#endif
