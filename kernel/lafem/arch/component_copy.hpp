// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_COPY_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_COPY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/half.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct ComponentCopy
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const int stride, const int block, const Index size)
        {
          BACKEND_SKELETON_VOID(value_cuda, value_generic, value_generic, r, x, stride, block, size)
        }

        template <typename DT_>
        static void value_to(const DT_ * const r, DT_ * x, const int stride, const int block, const Index size)
        {
          BACKEND_SKELETON_VOID(value_to_cuda, value_to_generic, value_to_generic, r, x, stride, block, size)
        }

        template <typename DT_>
        static void value_generic(DT_ * r, const DT_ * const x, const int stride, const int block, const Index size);

        template <typename DT_>
        static void value_to_generic(const DT_ * const r, DT_ * x, const int stride, const int block, const Index size);

        static void value_mkl(float * r, const float * const x, const int stride, const int block, const Index size);

        static void value_mkl(double * r, const double * const x, const int stride, const int block, const Index size);

        static void value_to_mkl(const float * const r, float * x, const int stride, const int block, const Index size);

        static void value_to_mkl(const double * const r, double * x, const int stride, const int block, const Index size);


        static void value_cuda(float * r, const float * const x, const int stride, const int block, const Index size);
        static void value_to_cuda(const float * const r, float * x, const int stride, const int block, const Index size);

        static void value_cuda(double * r, const double * const x, const int stride, const int block, const Index size);
        static void value_to_cuda(const double * const r, double * x, const int stride, const int block, const Index size);
#ifdef FEAT_HAVE_HALFMATH
        static void value_cuda(Half * r, const Half * const x, const int stride, const int block, const Index size);
        static void value_to_cuda(const Half * const r, Half * x, const int stride, const int block, const Index size);
#endif
      };

#ifdef FEAT_EICKT
      extern template void ComponentCopy::value_generic(float *, const float * const, const int, const int, const Index);
      extern template void ComponentCopy::value_generic(double *, const double * const, const int, const int, const Index);
      extern template void ComponentCopy::value_to_generic(const float * const, float *, const int, const int, const Index);
      extern template void ComponentCopy::value_to_generic(const double * const, double *, const int, const int, const Index);
#endif
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/component_copy_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_COMPONENT_COPY_HPP
