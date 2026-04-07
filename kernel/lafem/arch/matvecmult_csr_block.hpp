// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT::LAFEM::Arch
{
  /// r[i][k] <- y[i][k] + alpha * sum_{j=0..n-1} A[i][j] * x[j][k]
  struct MatVecMultCSRBlock
  {
    template <typename DT_, typename IT_, int block_size_>
    static void exec_generic_impl(Tiny::Vector<DT_,block_size_> * r, const DT_ alpha, const Tiny::Vector<DT_,block_size_> * const x, const Tiny::Vector<DT_,block_size_> * const y, const DT_ * const val,
      const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed);

    template <typename DT_, typename IT_, int block_size_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
      const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed);

    template <typename DT_, typename IT_, int block_size_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a,
      const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<Tiny::Vector<DT_, block_size_>> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<Tiny::Vector<DT_, block_size_>> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<Tiny::Vector<DT_, block_size_>> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl<DT_, IT_, block_size_>(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
      else
      {
        exec_generic_impl<DT_, IT_, block_size_>(r_view.get_w(), alpha, x_view.get_r(), (Tiny::Vector<DT_, block_size_>*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
    }

    template <typename DT_, typename IT_, int block_size_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a,
      const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl<DT_, IT_, block_size_>(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
      else
      {
        exec_cuda_impl<DT_, IT_, block_size_>(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, transposed);
      }
    }

    template <typename DT_, typename IT_, int block_size_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a,
      const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      XASSERT(r != x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_size_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MatVecMultCSRBlock::exec: MKL backend not available!");
        exec_generic<DT_, IT_, block_size_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
#else
        XABORTM("LAFEM::Arch::MatVecMultCSRBlock::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        if constexpr(block_size_ < 4)
          exec_cuda<DT_, IT_, block_size_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed);
        else
          XABORTM("LAFEM::Arch::MatVecMultCSRBlock::exec: CUDA backend not available for block_size_ > 3!");
#else
        XABORTM("LAFEM::Arch::MatVecMultCSRBlock::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultCSRBlock::exec: unknown backend!");
        break;
      }
    }
  }; // struct MatVecMultCSRBlock

#ifdef FEAT_EICKT
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint32_t, 1>(Tiny::Vector<float,1> *, const float, const Tiny::Vector<float,1> * const, const Tiny::Vector<float,1> * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint64_t, 1>(Tiny::Vector<float,1> *, const float, const Tiny::Vector<float,1> * const, const Tiny::Vector<float,1> * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint32_t, 2>(Tiny::Vector<float,2> *, const float, const Tiny::Vector<float,2> * const, const Tiny::Vector<float,2> * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint64_t, 2>(Tiny::Vector<float,2> *, const float, const Tiny::Vector<float,2> * const, const Tiny::Vector<float,2> * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint32_t, 3>(Tiny::Vector<float,3> *, const float, const Tiny::Vector<float,3> * const, const Tiny::Vector<float,3> * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint64_t, 3>(Tiny::Vector<float,3> *, const float, const Tiny::Vector<float,3> * const, const Tiny::Vector<float,3> * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint32_t, 4>(Tiny::Vector<float,4> *, const float, const Tiny::Vector<float,4> * const, const Tiny::Vector<float,4> * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<float, std::uint64_t, 4>(Tiny::Vector<float,4> *, const float, const Tiny::Vector<float,4> * const, const Tiny::Vector<float,4> * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint32_t, 1>(Tiny::Vector<double,1> *, const double, const Tiny::Vector<double,1> * const, const Tiny::Vector<double,1> * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint64_t, 1>(Tiny::Vector<double,1> *, const double, const Tiny::Vector<double,1> * const, const Tiny::Vector<double,1> * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint32_t, 2>(Tiny::Vector<double,2> *, const double, const Tiny::Vector<double,2> * const, const Tiny::Vector<double,2> * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint64_t, 2>(Tiny::Vector<double,2> *, const double, const Tiny::Vector<double,2> * const, const Tiny::Vector<double,2> * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint32_t, 3>(Tiny::Vector<double,3> *, const double, const Tiny::Vector<double,3> * const, const Tiny::Vector<double,3> * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint64_t, 3>(Tiny::Vector<double,3> *, const double, const Tiny::Vector<double,3> * const, const Tiny::Vector<double,3> * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint32_t, 4>(Tiny::Vector<double,4> *, const double, const Tiny::Vector<double,4> * const, const Tiny::Vector<double,4> * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
  extern template void MatVecMultCSRBlock::exec_generic_impl<double, std::uint64_t, 4>(Tiny::Vector<double,4> *, const double, const Tiny::Vector<double,4> * const, const Tiny::Vector<double,4> * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_CSR_BLOCK_HPP 1
#include <kernel/lafem/arch/matvecmult_csr_block_generic.hpp>
#endif
