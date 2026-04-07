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
  struct MatVecMultCSCRDense
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
      const IT_ * const col_idx, const IT_ * const row_ptr,  const IT_ * const row_idx, const Index, const Index, const Index, const Index, const bool);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Memory::Arbiter& row_idx, const Index rows, const Index cols, const Index nonzero_rows, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ri_view(row_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), ri_view.get_r(), rows, cols, nonzero_rows, nonzeros, transposed);
      }
      else
      {
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), ri_view.get_r(), rows, cols, nonzero_rows, nonzeros, transposed);
      }
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Memory::Arbiter& row_idx, const Index rows, const Index cols, const Index nonzero_rows, const Index nonzeros, const bool transposed)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(r, alpha, x, y, a, col_idx, row_ptr, row_idx, rows, cols, nonzero_rows, nonzeros, transposed);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MatVecMultCSCRDense::exec");
        exec_generic<DT_, IT_>(r, alpha, x, y, a, col_idx, row_ptr, row_idx, rows, cols, nonzero_rows, nonzeros, transposed);
#else
        XABORTM("LAFEM::Arch::MatVecMultCSCRDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
        XABORTM("LAFEM::Arch::MatVecMultCSCRDense::exec: CUDA backend not available!");
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultCSCRDense::exec: unknown backend!");
        break;
      }
    }
  }; // struct MatVecMultCSCRDense

#ifdef FEAT_EICKT
  extern template void MatVecMultCSCRDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSCRDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSCRDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSCRDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_CSCR_DENSE_HPP 1
#include <kernel/lafem/arch/matvecmult_cscr_dense_generic.hpp>
#endif
