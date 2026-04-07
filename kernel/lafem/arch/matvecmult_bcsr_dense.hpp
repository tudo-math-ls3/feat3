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
  /// r[i] <- y[i] + alpha * sum_{j=0..n-1} A[i][j] * x[j]
  struct MatVecMultBCSRDense
  {
    template <typename DT_, typename IT_, int bh_, int bw_>
    static void exec_generic_impl(Tiny::Vector<DT_, bh_> * r, const DT_ alpha, const Tiny::Vector<DT_, bw_> * const x,
      const Tiny::Vector<DT_, bh_> * const y, const Tiny::Matrix<DT_, bh_, bw_>* const val,
      const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index);

    template <typename DT_, typename IT_, int bh_, int bw_>
    static void exec_generic_transpose_impl(Tiny::Vector<DT_, bw_> * r, const DT_ alpha, const Tiny::Vector<DT_, bh_> * const x,
      const Tiny::Vector<DT_, bw_> * const y, const Tiny::Matrix<DT_, bh_, bw_>* const val,
      const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index);

    static void exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_size, const bool transposed);
    static void exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_size, const bool transposed);

    template <typename DT_, typename IT_>
    static void exec_mkl_impl(DT_ *, const DT_, const DT_ * const, const DT_ * const, const DT_ * const, const IT_ * const, const IT_ * const, const Index, const Index, const Index, const int, const bool)
    {
      XABORTM("LAFEM::Arch::MatVecMultBCSRDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const col_idx, const IT_ * const row_ptr, const Index rows, const Index cols, const Index nonzeros, const int block_height, const int block_width, const bool transposed);

    template <typename DT_, typename IT_, int block_height_, int block_width_, bool transposed_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::main));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::main));
      Memory::TypedView<Tiny::Matrix<DT_, block_height_, block_width_>> a_view(a.view(Memory::Location::main, Memory::Access::read));

      if constexpr (transposed_)
      {
        Memory::TypedView<Tiny::Vector<DT_, block_height_>> x_view(x.view(Memory::Location::main, Memory::Access::read));
        Memory::TypedView<Tiny::Vector<DT_, block_width_>> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
        if(y)
        {
          Memory::TypedView<Tiny::Vector<DT_, block_width_>> y_view(y.view(Memory::Location::main, Memory::Access::read));
          exec_generic_transpose_impl<DT_, IT_, block_height_, block_width_>(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
            ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros);
        }
        else
        {
          exec_generic_transpose_impl<DT_, IT_, block_height_, block_width_>(r_view.get_w(), alpha, x_view.get_r(), (Tiny::Vector<DT_, block_width_>*)nullptr, a_view.get_r(),
            ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros);
        }
      }
      else
      {
        Memory::TypedView<Tiny::Vector<DT_, block_width_>> x_view(x.view(Memory::Location::main, Memory::Access::read));
        Memory::TypedView<Tiny::Vector<DT_, block_height_>> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
        if(y)
        {
          Memory::TypedView<Tiny::Vector<DT_, block_height_>> y_view(y.view(Memory::Location::main, Memory::Access::read));
          exec_generic_impl<DT_, IT_, block_height_, block_width_>(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
            ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros);
        }
        else
        {
          exec_generic_impl<DT_, IT_, block_height_, block_width_>(r_view.get_w(), alpha, x_view.get_r(), (Tiny::Vector<DT_, block_height_>*)nullptr, a_view.get_r(),
            ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros);
        }
      }
    }

    template<typename DT_, int block_size_>
    static void exec_mkl(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<Index> rp_view(row_ptr.view(Memory::Location::main));
      Memory::TypedView<Index> ci_view(col_idx.view(Memory::Location::main));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, block_size_, transposed);
      }
      else
      {
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, block_size_, transposed);
      }
    }

    template<typename DT_, typename IT_, int block_height_, int block_width_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros, const bool transposed)
    {
      Memory::TypedView<IT_> rp_view(row_ptr.view(Memory::Location::cuda));
      Memory::TypedView<IT_> ci_view(col_idx.view(Memory::Location::cuda));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, block_height_, block_width_, transposed);
      }
      else
      {
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          ci_view.get_r(), rp_view.get_r(), rows, cols, nonzeros, block_height_, block_width_, transposed);
      }
    }

    template<typename DT_, typename IT_, int block_height_, int block_width_, bool transposed_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& col_idx, const Memory::Arbiter& row_ptr, const Index rows, const Index cols, const Index nonzeros)
    {
      XASSERT(r != x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_, block_height_, block_width_, transposed_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        if constexpr((block_height_ == block_width_) && (sizeof(IT_) == sizeof(Index)))
          exec_mkl<DT_, block_height_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed_);
        else
        {
          Backend::warn_missing("LAFEM::Arch::MatVecMultBCSRDense::exec: non-square blocks");
          exec_generic<DT_, IT_, block_height_, block_width_, transposed_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros);
        }
#else
        XABORTM("LAFEM::Arch::MatVecMultBCSRDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_, block_height_, block_width_>(r, alpha, x, y, a, col_idx, row_ptr, rows, cols, nonzeros, transposed_);
#else
        XABORTM("LAFEM::Arch::MatVecMultBCSRDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultBCSRDense::exec: unknown backend!");
        break;
      }
    }
  };

#ifdef FEAT_EICKT
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 1, 1>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 1, 1>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 2, 1>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 2, 1>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 3, 1>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 3, 1>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 4, 1>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 4, 1>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 1> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 1, 2>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 1, 2>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 2, 2>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 2, 2>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 3, 2>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 3, 2>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 4, 2>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 4, 2>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 2> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 1, 3>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 1, 3>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 2, 3>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 2, 3>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 3, 3>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 3, 3>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 4, 3>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 4, 3>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 3> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 1, 4>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 1, 4>(Tiny::Vector<float, 1>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 1> * const, const Tiny::Matrix<float, 1, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 2, 4>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 2, 4>(Tiny::Vector<float, 2>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 2> * const, const Tiny::Matrix<float, 2, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 3, 4>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 3, 4>(Tiny::Vector<float, 3>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 3> * const, const Tiny::Matrix<float, 3, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint32_t, 4, 4>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<float, std::uint64_t, 4, 4>(Tiny::Vector<float, 4>*, const float, const Tiny::Vector<float, 4> * const, const Tiny::Vector<float, 4> * const, const Tiny::Matrix<float, 4, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 1, 1>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 1, 1>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 2, 1>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 2, 1>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 3, 1>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 3, 1>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 4, 1>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 1>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 4, 1>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 1> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 1>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 1, 2>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 1, 2>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 2, 2>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 2, 2>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 3, 2>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 3, 2>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 4, 2>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 2>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 4, 2>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 2> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 2>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 1, 3>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 1, 3>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 2, 3>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 2, 3>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 3, 3>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 3, 3>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 4, 3>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 3>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 4, 3>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 3> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 3>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 1, 4>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 1, 4>(Tiny::Vector<double, 1>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 1> * const, const Tiny::Matrix<double, 1, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 2, 4>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 2, 4>(Tiny::Vector<double, 2>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 2> * const, const Tiny::Matrix<double, 2, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 3, 4>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 3, 4>(Tiny::Vector<double, 3>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 3> * const, const Tiny::Matrix<double, 3, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint32_t, 4, 4>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 4>* const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBCSRDense::exec_generic_impl<double, std::uint64_t, 4, 4>(Tiny::Vector<double, 4>*, const double, const Tiny::Vector<double, 4> * const, const Tiny::Vector<double, 4> * const, const Tiny::Matrix<double, 4, 4>* const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_BCSR_DENSE_HPP 1
#include <kernel/lafem/arch/matvecmult_bcsr_dense_generic.hpp>
#endif
