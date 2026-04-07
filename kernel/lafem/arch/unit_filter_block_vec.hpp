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
  struct UnitFilterBlockVec
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * v, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool zero, const bool ign_nans, const int block_size);

    template <typename DT_, typename IT_, int block_size_>
    static void exec_generic(Memory::Arbiter& v, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool zero, const bool ign_nans)
    {
      Memory::TypedView<DT_> v_view(v.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<IT_> i_view(f_idx.view(Memory::Location::main, Memory::Access::read));
      if(!zero || ign_nans)
      {
        Memory::TypedView<DT_> f_view(f_val.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(v_view.get_w(), f_view.get_r(), i_view.get_r(), f_nzes, zero, ign_nans, block_size_);
      }
      else
      {
        // f_val view is not required in this case
        exec_generic_impl(v_view.get_w(), (DT_*)nullptr, i_view.get_r(), f_nzes, zero, ign_nans, block_size_);
      }
    }

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * v, const DT_ * const f_val, const IT_ * const f_idx, const Index f_nzes, const bool zero, const bool ign_nans, const int block_size);

    template <typename DT_, typename IT_, int block_size_>
    static void exec_cuda(Memory::Arbiter& v, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool zero, const bool ign_nans)
    {
      Memory::TypedView<DT_> v_view(v.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<IT_> i_view(f_idx.view(Memory::Location::cuda, Memory::Access::read));
      if(!zero || ign_nans)
      {
        Memory::TypedView<DT_> f_view(f_val.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(v_view.get_w(), f_view.get_r(), i_view.get_r(), f_nzes, zero, ign_nans, block_size_);
      }
      else
      {
        // f_val view is not required in this case
        exec_cuda_impl(v_view.get_w(), (DT_*)nullptr, i_view.get_r(), f_nzes, zero, ign_nans, block_size_);
      }
    }

    template <typename DT_, typename IT_, int block_size_>
    static void exec(Memory::Arbiter& v, const Memory::Arbiter& f_val, const Memory::Arbiter& f_idx, const Index f_nzes, const bool zero, const bool ign_nans)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_size_>(v, f_val, f_idx, f_nzes, zero, ign_nans);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::UnitFilterBlockVec::exec");
        exec_generic<DT_, IT_, block_size_>(v, f_val, f_idx, f_nzes, zero, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockVec::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, block_size_>(v, f_val, f_idx, f_nzes, zero, ign_nans);
#else
        XABORTM("LAFEM::Arch::UnitFilterBlockVec::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::UnitFilterBlockVec::exec: unknown backend!");
        break;
      }
    }
  }; // struct UnitFilterBlockVec

#ifdef FEAT_EICKT
  extern template void UnitFilterBlockVec::exec_generic_impl(float *, const float * const, const std::uint64_t * const, const Index, const bool, const bool, const int);
  extern template void UnitFilterBlockVec::exec_generic_impl(double *, const double * const, const std::uint64_t * const, const Index, const bool, const bool, const int);
  extern template void UnitFilterBlockVec::exec_generic_impl(float *, const float * const, const std::uint32_t * const, const Index, const bool, const bool, const int);
  extern template void UnitFilterBlockVec::exec_generic_impl(double *, const double * const, const std::uint32_t * const, const Index, const bool, const bool, const int);
#endif
} // namespace FEAT::LAFEM::Arch
/// \endcond

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_BLOCK_VEC_HPP 1
#include <kernel/lafem/arch/unit_filter_block_vec_generic.hpp>
#endif
