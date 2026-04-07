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
  struct MatMatMultSparseDense
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index nonzeros,
      const DT_ * const y, const DT_* const z, const Index rows,  const Index cols, const Index inner);

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index nonzeros,
      const DT_ * const y, const DT_* const z, const Index rows,  const Index cols, const Index inner);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index nonzeros,
      const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows,  const Index cols, const Index inner)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(z)
      {
        Memory::TypedView<DT_> z_view(z.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, a_view.get_r(), ci_view.get_r(), rp_view.get_r(),
          nonzeros, y_view.get_r(), z_view.get_r(), rows, cols, inner);
      }
      else
      {
        exec_generic_impl(r_view.get_w(), alpha, a_view.get_r(), ci_view.get_r(), rp_view.get_r(),
          nonzeros, y_view.get_r(), (DT_*)nullptr, rows, cols, inner);
      }
    }

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index nonzeros,
      const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows,  const Index cols, const Index inner)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      if(z)
      {
        Memory::TypedView<DT_> z_view(z.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, a_view.get_r(), ci_view.get_r(), rp_view.get_r(),
          nonzeros, y_view.get_r(), z_view.get_r(), rows, cols, inner);
      }
      else
      {
        exec_cuda_impl(r_view.get_w(), alpha, a_view.get_r(), ci_view.get_r(), rp_view.get_r(),
          nonzeros, y_view.get_r(), (DT_*)nullptr, rows, cols, inner);
      }
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& val, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index nonzeros,
      const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows,  const Index cols, const Index inner)
    {
      XASSERT(r != y);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(r, alpha, val, col_idx, row_ptr, nonzeros, y, z, rows, cols, inner);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MatMatMultSparseDense::exec");
        exec_generic<DT_, IT_>(r, alpha, val, col_idx, row_ptr, nonzeros, y, z, rows, cols, inner);
#else
        XABORTM("LAFEM::Arch::MatMatMultSparseDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(r, alpha, val, col_idx, row_ptr, nonzeros, y, z, rows, cols, inner);
#else
        XABORTM("LAFEM::Arch::MatMatMultSparseDense::exec: CUDA backend not available!");
#endif
        break;
      }
    }
  }; // struct MatMatMultSparseDense

#ifdef FEAT_EICKT
  extern template void MatMatMultSparseDense::exec_generic_impl(float *, const float, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const float * const, const float * const, const Index, const Index, const Index);
  extern template void MatMatMultSparseDense::exec_generic_impl(double *, const double, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const double *const , const double *const , const Index, const Index, const Index);

  extern template void MatMatMultSparseDense::exec_generic_impl(float *, const float, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const float * const, const float * const, const Index, const Index, const Index);
  extern template void MatMatMultSparseDense::exec_generic_impl(double *, const double, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const double * const, const double *const , const Index, const Index, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATMATMULT_SPARSE_DENSE_HPP 1
#include <kernel/lafem/arch/matmatmult_sparse_dense_generic.hpp>
#endif
