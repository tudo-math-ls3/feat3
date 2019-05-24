// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
#define KERNEL_LAFEM_ARCH_AXPY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/archs.hpp>

#include <typeinfo>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct Axpy;

      template <>
      struct Axpy<Mem::Main>
      {
        template <typename DT_>
        static void dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
        {
          dv_generic(r, a, x, y, size);
        }

#ifdef FEAT_HAVE_MKL
        static void dv(float * r, const float a, const float * const x, const float * const y, const Index size)
        {
        /// \compilerhack icc (in combination with mkl)  crashes kernel/solver/optimiser-test when calculating axpy on very small vectors (size=2)
#ifdef FEAT_COMPILER_INTEL
          if (size < 17)
            dv_generic(r, a, x, y, size);
          else
#endif
            dv_mkl(r, a, x, y, size);
        }

        static void dv(double * r, const double a, const double * const x, const double * const y, const Index size)
        {
        /// \compilerhack icc (in combination with mkl)  crashes kernel/solver/optimiser-test when calculating axpy on very small vectors (size=2)
#ifdef FEAT_COMPILER_INTEL
          if (size < 17)
            dv_generic(r, a, x, y, size);
          else
#endif
            dv_mkl(r, a, x, y, size);
        }
#endif

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void dv(__float128 * r, const __float128 a, const __float128 * const x, const __float128 * const y, const Index size)
        {
          dv_generic(r, a, x, y, size);
        }
#endif

        template <typename DT_>
        static void dv_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);

        static void dv_mkl(float * r, const float a, const float * const x, const float * const y, const Index size);
        static void dv_mkl(double * r, const double a, const double * const x, const double * const y, const Index size);
      };

#ifdef FEAT_EICKT
      extern template void Axpy<Mem::Main>::dv_generic(float *, const float, const float * const, const float * const, const Index);
      extern template void Axpy<Mem::Main>::dv_generic(double *, const double, const double * const, const double * const, const Index);
#endif

      template <>
      struct Axpy<Mem::CUDA>
      {
        template <typename DT_>
        static void dv(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/axpy_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_AXPY_HPP
