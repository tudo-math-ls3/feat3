// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_AXPY_HPP
#define KERNEL_LAFEM_ARCH_AXPY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Axpy
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size)
        {
          value_generic(r, a, x, y, size);
        }

#ifdef FEAT_HAVE_HALFMATH
        static void value(Half * r, const Half a, const Half * const x, const Half * const y, const Index size)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_generic, value_generic, r, a, x, y, size)
        }
#endif

        static void value(float * r, const float a, const float * const x, const float * const y, const Index size)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_mkl, value_generic, r, a, x, y, size)
        }

        static void value(double * r, const double a, const double * const x, const double * const y, const Index size)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_mkl, value_generic, r, a, x, y, size)
        }

        template <typename DT_>
        static void value_generic(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);

        static void value_mkl(float * r, const float a, const float * const x, const float * const y, const Index size);
        static void value_mkl(double * r, const double a, const double * const x, const double * const y, const Index size);

        template <typename DT_>
        static void value_cuda(DT_ * r, const DT_ a, const DT_ * const x, const DT_ * const y, const Index size);
      };

#ifdef FEAT_EICKT
      extern template void Axpy::value_generic(float *, const float, const float * const, const float * const, const Index);
      extern template void Axpy::value_generic(double *, const double, const double * const, const double * const, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/axpy_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_AXPY_HPP
