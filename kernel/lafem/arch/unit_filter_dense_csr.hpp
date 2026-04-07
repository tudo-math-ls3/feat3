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
  struct UnitFilterDenseCSR
  {
    template<typename DT_, typename IT_>
    static void exec_generic_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const f_idx, const Index f_nzes, const bool unit, const int block_width);

    template<typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const int block_width)
    {
      Memory::TypedView<DT_> a_view(mat.view(Memory::Location::main, Memory::Access::read_write));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::main, Memory::Access::read));
      exec_generic_impl(a_view.get_w(), rp_view.get_r(), ci_view.get_r(), fi_view.get_r(), f_nzes, unit, block_width);
    }

    template<typename DT_, typename IT_>
    static void exec_cuda_impl(DT_* mat, const IT_* const row_ptr, const IT_* const col_idx, const IT_ * const f_idx, const Index f_nzes, const bool unit, const int block_width);

    template<typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const int block_width)
    {
      Memory::TypedView<DT_> a_view(mat.view(Memory::Location::cuda, Memory::Access::read_write));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> fi_view(f_idx.view(Memory::Location::cuda, Memory::Access::read));
      exec_cuda_impl(a_view.get_w(), rp_view.get_r(), ci_view.get_r(), fi_view.get_r(), f_nzes, unit, block_width);
    }

    template<typename DT_, typename IT_>
    static void exec(Memory::Arbiter& mat, const Memory::Arbiter& row_ptr, const Memory::Arbiter& col_idx, const Memory::Arbiter& f_idx, const Index f_nzes, const bool unit, const int block_width)
    {
      XASSERT(!unit || (block_width == 1));

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(mat, row_ptr, col_idx, f_idx, f_nzes, unit, block_width);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::UnitFilterDenseCSR::exec");
        exec_generic<DT_, IT_>(mat, row_ptr, col_idx, f_idx, f_nzes, unit, block_width);
#else
        XABORTM("LAFEM::Arch::UnitFilterDenseCSR::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(mat, row_ptr, col_idx, f_idx, f_nzes, unit, block_width);
#else
        XABORTM("LAFEM::Arch::UnitFilterDenseCSR::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::UnitFilterDenseCSR::exec: unknown backend!");
        break;
      }
    }
  }; // struct UnitFilterDenseCSR

#ifdef FEAT_EICKT
  extern template void UnitFilterDenseCSR::exec_generic_impl(float* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const int);
  extern template void UnitFilterDenseCSR::exec_generic_impl(float* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const int);
  extern template void UnitFilterDenseCSR::exec_generic_impl(double* mat, const std::uint32_t* const row_ptr, const std::uint32_t* const col_idx, const std::uint32_t * const f_idx, const Index f_nzes, const bool unit, const int);
  extern template void UnitFilterDenseCSR::exec_generic_impl(double* mat, const std::uint64_t* const row_ptr, const std::uint64_t* const col_idx, const std::uint64_t * const f_idx, const Index f_nzes, const bool unit, const int);
#endif
} // namespace FEAT::LAFEM::Arch
/// \endcond

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_UNIT_FILTER_DENSE_CSR_HPP 1
#include <kernel/lafem/arch/unit_filter_dense_csr_generic.hpp>
#endif
