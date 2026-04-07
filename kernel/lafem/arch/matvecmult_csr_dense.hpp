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
  /// r[i] <- y[i] + alpha * sum_{j=0..n-1} A[i][j] * x[j]
  struct MatVecMultCSRDense
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
      const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index, const Index, const bool);

    static void exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const bool);
    static void exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const bool);

    template <typename DT_, typename IT_>
    static void exec_mkl_impl(DT_ *, const DT_, const DT_ * const, const DT_ * const, const DT_ * const, const IT_ * const, const IT_ * const, const Index, const Index, const Index, const bool)
    {
      XABORTM("LAFEM::Arch::MatVecMultCSRDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr,  const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
       }
      else // no y
      {
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
    }

    template <typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      // we only expect MKL to work with indices of the same size as FEAT::Index
      Memory::TypedView<Index> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<Index> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
      else // no y
      {
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
    }

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
      else
      {
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      XASSERT(r != x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        if constexpr(sizeof(IT_) == sizeof(Index))
          exec_mkl<DT_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
        else
        {
          Backend::warn_missing("LAFEM::Arch::MatVecMultCSRDense::exec: MKL backend not available for IT_ != Index");
          exec_generic<DT_, IT_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
        }
#else
        XABORTM("LAFEM::Arch::MatVecMultCSRDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
#else
        XABORTM("LAFEM::Arch::MatVecMultCSRDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultCSRDense::exec: unknown backend!");
        break;
      }
    }
  }; // struct MatVecMultCSRDense

#ifdef FEAT_EICKT
  extern template void MatVecMultCSRDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_CSR_DENSE_HPP 1
#include <kernel/lafem/arch/matvecmult_csr_dense_generic.hpp>
#endif
