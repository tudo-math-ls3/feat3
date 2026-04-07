// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_arbiter.hpp>

/// \cond internal
namespace FEAT::LAFEM::Arch
{
  struct SlipFilterBlock
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * v, const DT_ * const f_nu, const IT_ * const f_idx, const Index f_nzes, const int bs);

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * v, const DT_ * const f_nu, const IT_ * const f_idx, const Index f_nzes, const int bs);

    template <typename DT_, typename IT_, int block_size_>
    static void exec_generic(Memory::Arbiter& v, const Memory::Arbiter& f_nu, const Memory::Arbiter& f_idx, const Index f_nzes)
    {
      Memory::TypedView<DT_> v_view(v.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<DT_> n_view(f_nu.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> i_view(f_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(v_view.get_w(), n_view.get_r(), i_view.get_r(), f_nzes, block_size_);
    }

    template <typename DT_, typename IT_, int block_size_>
    static void exec_cuda(Memory::Arbiter& v, const Memory::Arbiter& f_nu, const Memory::Arbiter& f_idx, const Index f_nzes)
    {
      Memory::TypedView<DT_> v_view(v.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<DT_> n_view(f_nu.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> i_view(f_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(v_view.get_w(), n_view.get_r(), i_view.get_r(), f_nzes, block_size_);
    }

    template <typename DT_, typename IT_, int block_size_>
    static void exec(Memory::Arbiter& v, const Memory::Arbiter& f_nu, const Memory::Arbiter& f_idx, const Index f_nzes)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_size_>(v, f_nu, f_idx, f_nzes);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::SlipFilterBlock::exec");
        exec_generic<DT_, IT_, block_size_>(v, f_nu, f_idx, f_nzes);
#else
        XABORTM("LAFEM::Arch::SlipFilterBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, block_size_>(v, f_nu, f_idx, f_nzes);
#else
        XABORTM("LAFEM::Arch::SlipFilterBlock::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::SlipFilterBlock::exec: unknown backend!");
        break;
      }
    }
  }; // SlipFilterBlock

#ifdef FEAT_EICKT
  extern template void SlipFilterBlock::exec_generic_impl<float, std::uint64_t>(float * v, const float * const f_nu, const std::uint64_t * const f_idx, const Index f_nzes, const int bs);
  extern template void SlipFilterBlock::exec_generic_impl<double, std::uint64_t>(double * v, const double * const f_nu, const std::uint64_t * const f_idx, const Index f_nzes, const int bs);
  extern template void SlipFilterBlock::exec_generic_impl<float, std::uint32_t>(float * v, const float * const f_nu, const std::uint32_t * const f_idx, const Index f_nzes, const int bs);
  extern template void SlipFilterBlock::exec_generic_impl<double, std::uint32_t>(double * v, const double * const f_nu, const std::uint32_t * const f_idx, const Index f_nzes, const int bs);
#endif
} // namespace FEAT::LAFEM::Arch

/// \endcond
#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_SLIP_FILTER_BLOCK_HPP 1
#include <kernel/lafem/arch/slip_filter_block_generic.hpp>
#endif
