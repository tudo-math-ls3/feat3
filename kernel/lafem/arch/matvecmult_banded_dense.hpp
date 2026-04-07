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
  struct MatVecMultBandedDense
  {
    template <typename DT_, typename IT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val,
      const IT_ * const offsets, const Index bands, const Index rows, const Index cols);

    template <typename DT_, typename IT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const val, const IT_ * const offsets, const Index bands, const Index rows, const Index cols);

    template <typename DT_, typename IT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& offsets, const Index bands, const Index rows, const Index cols)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<IT_> o_view(offsets.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          o_view.get_r(), bands, rows, cols);
      }
      else
      {
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          o_view.get_r(), bands, rows, cols);
      }
    }

    template <typename DT_, typename IT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& offsets, const Index bands, const Index rows, const Index cols)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<IT_> o_view(offsets.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> a_view(a.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), a_view.get_r(),
          o_view.get_r(), bands, rows, cols);
      }
      else
      {
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), (DT_*)nullptr, a_view.get_r(),
          o_view.get_r(), bands, rows, cols);
      }
    }

    template <typename DT_, typename IT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& a, const Memory::Arbiter& offsets, const Index bands, const Index rows, const Index cols)
    {
      XASSERT(r != x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_, IT_>(r, alpha, x, y, a, offsets, bands, rows, cols);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        Backend::warn_missing("LAFEM::Arch::MatVecMultBandedDense::exec");
        exec_generic<DT_, IT_>(r, alpha, x, y, a, offsets, bands, rows, cols);
#else
        XABORTM("LAFEM::Arch::MatVecMultBandedDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, IT_>(r, alpha, x, y, a, offsets, bands, rows, cols);
#else
        XABORTM("LAFEM::Arch::MatVecMultBandedDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultBandedDense::exec: unknown backend!");
        break;
      }
    }
  }; // struct MatVecMultBandedDense

#ifdef FEAT_EICKT
  extern template void MatVecMultBandedDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBandedDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const std::uint32_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBandedDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint64_t * const, const Index, const Index, const Index);
  extern template void MatVecMultBandedDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const std::uint32_t * const, const Index, const Index, const Index);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_BANDED_DENSE_HPP 1
#include <kernel/lafem/arch/matvecmult_banded_dense_generic.hpp>
#endif
