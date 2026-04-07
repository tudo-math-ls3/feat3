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
  struct MirrorDenseGather
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

    template<typename DT_, typename IT_>
    static void exec_cuda_impl(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

    template<typename DT_, typename IT_>
    static void exec_generic(const Index boff, const Index nidx, const Memory::Arbiter& idx, Memory::Arbiter& buf, const Memory::Arbiter& vec)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::main, Memory::Access::write));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::main, Memory::Access::read));
      return exec_generic_impl(boff, nidx, i_view.get_r(), b_view.get_w(), x_view.get_r());
    }

    template<typename DT_, typename IT_>
    static void exec_cuda(const Index boff, const Index nidx, const Memory::Arbiter& idx, Memory::Arbiter& buf, const Memory::Arbiter& vec)
    {
      Memory::TypedView<DT_> b_view(buf.view(Memory::Location::cuda, Memory::Access::write));
      Memory::TypedView<DT_> x_view(vec.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> i_view(idx.view(Memory::Location::cuda, Memory::Access::read));
      return exec_cuda_impl(boff, nidx, i_view.get_r(), b_view.get_w(), x_view.get_r());
    }

    template<typename DT_, typename IT_>
    static void exec(const Index boff, const Index nidx, const Memory::Arbiter& idx, Memory::Arbiter& buf, const Memory::Arbiter& vec)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(boff, nidx, idx, buf, vec);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MirrorDenseGather::exec");
        exec_generic<DT_, IT_>(boff, nidx, idx, buf, vec);
#else
        XABORTM("LAFEM::Arch::MirrorDenseGather::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(boff, nidx, idx, buf, vec);
#else
        XABORTM("LAFEM::Arch::MirrorDenseGather::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MirrorDenseGather::exec: unknown backend!");
        break;
      }
    }
  }; // struct Mirror
} // namespace FEAT::LAFEM::Arch

#ifndef __CUDACC__
#define KERNEL_LAFEM_ARCH_MIRROR_DENSE_GATHER_HPP 1
#include <kernel/lafem/arch/mirror_dense_gather_generic.hpp>
#endif
