// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP 1

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
      struct DotProduct
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size)
        {
          return value_generic(x, y, size);
        }

        template <typename ValueType_>
        static ValueType_ value_blocked(const ValueType_ * const x, const ValueType_ * const y, const Index size)
        {
          return value_blocked_generic(x, y, size);
        }

#ifdef FEAT_HAVE_HALFMATH
        static Half value(const Half * const x, const Half * const y, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_generic, value_generic, x, y, size)
        }
#endif

        static float value(const float * const x, const float * const y, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, y, size)
        }

        static double value(const double * const x, const double * const y, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, y, size)
        }

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename ValueType_>
        static ValueType_ value_blocked_generic(const ValueType_ * const x, const ValueType_ * const y, const Index size);

        static float value_mkl(const float * const x, const float * const y, const Index size);
        static double value_mkl(const double * const x, const double * const y, const Index size);

        template <typename DT_>
        static DT_ value_cuda(const DT_ * const x, const DT_ * const y, const Index size);
      };

#ifdef FEAT_EICKT
      extern template float DotProduct::value_generic(const float * const, const float * const, const Index);
      extern template double DotProduct::value_generic(const double * const, const double * const, const Index);
#endif

      struct TripleDotProduct
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
        {
          return value_generic(x, y, z, size);
        }

        template <typename ValueType_>
        static ValueType_ value_blocked(const ValueType_ * const x, const ValueType_ * const y, const ValueType_ * const z, const Index size)
        {
          return value_blocked_generic(x, y, z, size);
        }

#ifdef FEAT_HAVE_HALFMATH
        static Half value(const Half * const x, const Half * const y, const Half * const z, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_generic, value_generic, x, y, z, size)
        }
#endif

        static float value(const float * const x, const float * const y, const float * const z, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, y, z, size)
        }

        static double value(const double * const x, const double * const y, const double * const z, const Index size)
        {
          BACKEND_SKELETON_RETURN(value_cuda, value_mkl, value_generic, x, y, z, size)
        }

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);

        template <typename ValueType_>
        static ValueType_ value_blocked_generic(const ValueType_ * const x, const ValueType_ * const y, const ValueType_ * const z, const Index size);

        static float value_mkl(const float * const x, const float * const y, const float * const z, const Index size);
        static double value_mkl(const double * const x, const double * const y, const double * const z, const Index size);

        template <typename DT_>
        static DT_ value_cuda(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

#ifdef FEAT_EICKT
      extern template float TripleDotProduct::value_generic(const float * const, const float * const, const float * const, const Index);
      extern template double TripleDotProduct::value_generic(const double * const, const double * const, const double * const, const Index);
#endif
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/dot_product_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
