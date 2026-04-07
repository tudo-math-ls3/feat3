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
  struct MatVecMultDenseDense
  {
    template <typename DT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index cols, bool transpose);

    static void exec_mkl_impl(float * r, const float alpha, const float * const y, const float * const val, const float * const x, const Index rows, const Index cols, bool transpose);
    static void exec_mkl_impl(double * r, const double alpha, const double * const y, const double * const val, const double * const x, const Index rows, const Index cols, bool transpose);

    template <typename DT_>
    static void exec_mkl_impl(DT_ *, const DT_, const DT_ * const, const DT_ * const, const DT_ * const, const Index, const Index, bool)
    {
      XABORTM("LAFEM::Arch::MatVecMultDenseDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha, const DT_ * const y, const DT_ * const val, const DT_ * const x, const Index rows, const Index cols, bool transpose);

    template <typename DT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& y, const Memory::Arbiter& val, const Memory::Arbiter& x, const Index rows, const Index cols, bool transposed)
    {
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, y_view.get_r(), a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
      else
      {
        exec_generic_impl(r_view.get_w(), alpha, (DT_*)nullptr, a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
    }

    template <typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& y, const Memory::Arbiter& val, const Memory::Arbiter& x, const Index rows, const Index cols, bool transposed)
    {
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
        exec_mkl_impl(r_view.get_w(), alpha, y_view.get_r(), a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
      else
      {
        exec_mkl_impl(r_view.get_w(), alpha, (DT_*)nullptr, a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
    }

    template <typename DT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& y, const Memory::Arbiter& val, const Memory::Arbiter& x, const Index rows, const Index cols, bool transposed)
    {
      Memory::TypedView<DT_> a_view(val.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      if(y)
      {
        Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, y_view.get_r(), a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
      else
      {
        exec_cuda_impl(r_view.get_w(), alpha, (DT_*)nullptr, a_view.get_r(), x_view.get_r(), rows, cols, transposed);
      }
    }

    template <typename DT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& y, const Memory::Arbiter& val, const Memory::Arbiter& x, const Index rows, const Index cols, bool transposed)
    {
      XASSERT(r != x);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic(r, alpha, y, val, x, rows, cols, transposed);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        exec_mkl(r, alpha, y, val, x, rows, cols, transposed);
#else
        XABORTM("LAFEM::Arch::MatVecMultDenseDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda(r, alpha, y, val, x, rows, cols, transposed);
#else
        XABORTM("LAFEM::Arch::MatVecMultDenseDense::exec: CUDA backend not available!");
#endif
        break;

      default:
        XABORTM("LAFEM::Arch::MatVecMultDenseDense::exec: unknown backend!");
        break;
      }
    }
  }; // struct MatVecMultDenseDense

#ifdef FEAT_EICKT
  extern template void MatVecMultDenseDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const Index, const Index, bool);
  extern template void MatVecMultDenseDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const Index, const Index, bool);
#endif

} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATVECMULT_DENSE_DENSE_HPP 1
#include <kernel/lafem/arch/matvecmult_dense_dense_generic.hpp>
#endif
