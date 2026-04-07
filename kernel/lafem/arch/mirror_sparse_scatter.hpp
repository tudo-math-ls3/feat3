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
  struct MirrorSparseScatter
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha);

    template<typename DT_, typename IT_>
    static void exec_generic(const Index boff, const Index nidx, const Memory::Arbiter& idx, const Memory::Arbiter& buf, const Index nvec, Memory::Arbiter& vec, const Memory::Arbiter& vi, const DT_ alpha)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> vi_view(vi.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(boff, nidx, i_view.get_r(), b_view.get_r(), nvec, x_view.get_w(), vi_view.get_r(), alpha);
    }

    template<typename DT_, typename IT_>
    static void exec(const Index boff, const Index nidx, const Memory::Arbiter& idx, const Memory::Arbiter& buf, const Index nvec, Memory::Arbiter& vec, const Memory::Arbiter& vi, const DT_ alpha)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        return exec_generic<DT_, IT_>(boff, nidx, idx, buf, nvec, vec, vi, alpha);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MirrorSparseScatter::exec");
        return exec_generic<DT_, IT_>(boff, nidx, idx, buf, nvec, vec, vi, alpha);
#else
        XABORTM("LAFEM::Arch::MirrorSparseScatter::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::MirrorSparseScatter::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::MirrorSparseScatter::exec: unknown backend!");
        break;
      }
    }
  }; // struct MirrorSparseScatter
} // namespace FEAT::LAFEM::Arch

#ifndef __CUDACC__
#define KERNEL_LAFEM_ARCH_MIRROR_SPARSE_SCATTER_HPP 1
#include <kernel/lafem/arch/mirror_sparse_scatter_generic.hpp>
#endif
