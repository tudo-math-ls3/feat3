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
  struct MirrorSparseBlockGather
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx);

    template<typename DT_, typename IT_, int block_size_>
    static void exec_generic(const Index boff, const Index nidx, const Memory::Arbiter& idx, Memory::Arbiter& buf, Index nvec, const Memory::Arbiter& vec, const Memory::Arbiter& vi)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> vi_view(vi.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(Index(block_size_), boff, nidx, i_view.get_r(), b_view.get_w(), nvec, x_view.get_r(), vi_view.get_r());
    }

    template<typename DT_, typename IT_, int block_size_>
    static void exec(const Index boff, const Index nidx, const Memory::Arbiter& idx, Memory::Arbiter& buf, Index nvec, const Memory::Arbiter& vec, const Memory::Arbiter& vi)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_size_>(boff, nidx, idx, buf, nvec, vec, vi);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MirrorSparseBlockGather::exec");
        exec_generic<DT_, IT_, block_size_>(boff, nidx, idx, buf, nvec, vec, vi);
#else
        XABORTM("LAFEM::Arch::MirrorSparseBlockGather::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::MirrorSparseBlockGather::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::MirrorSparseBlockGather::exec: unknown backend!");
        break;
      }
    }
  }; // struct MirrorSparseBlockGather
} // namespace FEAT::LAFEM::Arch

#ifndef __CUDACC__
#define KERNEL_LAFEM_ARCH_MIRROR_SPARSE_BLOCK_GATHER_HPP 1
#include <kernel/lafem/arch/mirror_sparse_block_gather_generic.hpp>
#endif
