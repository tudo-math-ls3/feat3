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
  struct MatMatMultDenseDense
  {
    template <typename DT_>
    static void exec_generic_impl(DT_ * r, const DT_ alpha, const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index cols, const Index inner);

    static void exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const z, const Index rows, const Index cols, const Index inner);
    static void exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const z, const Index rows, const Index cols, const Index inner);

    template <typename DT_>
    static void exec_mkl_impl(DT_ *, const DT_,  const DT_ * const, const DT_ * const, const DT_ * const, const Index, const Index, const Index)
    {
      XABORTM("LAFEM::Arch::MatMatMultDenseDense::exec_mkl_impl: MKL backend not available!");
    }

    template <typename DT_>
    static void exec_cuda_impl(DT_ * r, const DT_ alpha,  const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index rows, const Index cols, const Index inner);

    template <typename DT_>
    static void exec_generic(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows, const Index cols, const Index inner)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      if(z)
      {
        Memory::TypedView<DT_> z_view(z.view(Memory::Location::main, Memory::Access::read));
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), z_view.get_r(), rows, cols, inner);
      }
      else
        exec_generic_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), (DT_*)nullptr, rows, cols, inner);
    }

    template <typename DT_>
    static void exec_mkl(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows, const Index cols, const Index inner)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::main, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::main, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::main, Memory::Access::read));
      if(z)
      {
        Memory::TypedView<DT_> z_view(z.view(Memory::Location::main, Memory::Access::read));
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), z_view.get_r(), rows, cols, inner);
      }
      else
        exec_mkl_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), (DT_*)nullptr, rows, cols, inner);
    }

    template <typename DT_>
    static void exec_cuda(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows, const Index cols, const Index inner)
    {
      Memory::TypedView<DT_> r_view(r.view(Memory::Location::cuda, Memory::Access::write | Memory::Access::overlap));
      Memory::TypedView<DT_> x_view(x.view(Memory::Location::cuda, Memory::Access::read));
      Memory::TypedView<DT_> y_view(y.view(Memory::Location::cuda, Memory::Access::read));
      if(z)
      {
        Memory::TypedView<DT_> z_view(z.view(Memory::Location::cuda, Memory::Access::read));
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), z_view.get_r(), rows, cols, inner);
      }
      else
        exec_cuda_impl(r_view.get_w(), alpha, x_view.get_r(), y_view.get_r(), (DT_*)nullptr, rows, cols, inner);
    }

    template <typename DT_>
    static void exec(Memory::Arbiter& r, const DT_ alpha, const Memory::Arbiter& x, const Memory::Arbiter& y, const Memory::Arbiter& z, const Index rows, const Index cols, const Index inner)
    {
      XASSERT(r != x);
      XASSERT(r != y);

      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::generic:
        exec_generic<DT_>(r, alpha, x, y, z, rows, cols, inner);
        break;

      case PreferredBackend::mkl:
#ifdef FEAT_HAVE_MKL
        exec_mkl<DT_>(r, alpha, x, y, z, rows, cols, inner);
#else
        XABORTM("LAFEM::Arch::MatMatMultDenseDense::exec: MKL backend not available!");
#endif
        break;

      case PreferredBackend::cuda:
#ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_>(r, alpha, x, y, z, rows, cols, inner);
#else
        XABORTM("LAFEM::Arch::MatMatMultDenseDense::exec: CUDA backend not available!");
#endif
        break;
      }
    }
  }; // struct MatMatMultDenseDense

#ifdef FEAT_EICKT
  extern template void MatMatMultDenseDense::exec_generic_impl(float *, const float, const float * const, const float * const, const float * const, const Index, const Index, const Index);
  extern template void MatMatMultDenseDense::exec_generic_impl(double *, const double, const double * const, const double * const, const double * const, const Index, const Index, const Index);
#endif
} // namespace FEAT::LAFEM::Arch

#ifndef  __CUDACC__
#define KERNEL_LAFEM_ARCH_MATMATMULT_DENSE_DENSE_HPP 1
#include <kernel/lafem/arch/matmatmult_dense_dense_generic.hpp>
#endif
