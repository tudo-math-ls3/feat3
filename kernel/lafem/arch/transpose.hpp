// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_TRANSPOSE_HPP
#define KERNEL_LAFEM_ARCH_TRANSPOSE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Transpose
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const Index rows_x, const Index columns_x)
        {
          value_generic(r, x, rows_x, columns_x);
        }

        static void value(float * r, const float * const x, const Index rows_x, const Index columns_x)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_mkl, value_generic, r, x, rows_x, columns_x)
        }

        static void value(double * r, const double * const x, const Index rows_x, const Index columns_x)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_mkl, value_generic, r, x, rows_x, columns_x)
        }

        template <typename DT_>
        static void value_generic(DT_ * r, const DT_ * const x, const Index rows_x, const Index columns_x);

        static void value_mkl(float * r, const float * const x, const Index rows_x, const Index columns_x);
        static void value_mkl(double * r, const double * const x, const Index rows_x, const Index columns_x);

        static void value_cuda(float * r, const float * const x, const Index rows_x, const Index columns_x);
        static void value_cuda(double * r, const double * const x, const Index rows_x, const Index columns_x);
      };

#ifdef FEAT_EICKT
      extern template void Transpose::value_generic(float *, const float * const, const Index, const Index);
      extern template void Transpose::value_generic(double *, const double * const, const Index, const Index);
#endif

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/transpose_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_TRANSPOSE_HPP
