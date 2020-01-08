// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>



namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct DotProduct;

      template <>
      struct DotProduct<Mem::Main>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size)
        {
          return value_generic(x, y, size);
        }

#ifdef FEAT_HAVE_MKL
        static float value(const float * const x, const float * const y, const Index size)
        {
          return value_mkl(x, y, size);
        }

        static double value(const double * const x, const double * const y, const Index size)
        {
          return value_mkl(x, y, size);
        }
#endif // FEAT_HAVE_MKL

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
        static __float128 value(const __float128 * const x, const __float128 * const y, const Index size)
        {
          return value_generic(x, y, size);
        }
#endif

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const DT_ * const y, const Index size);

        static float value_mkl(const float * const x, const float * const y, const Index size);
        static double value_mkl(const double * const x, const double * const y, const Index size);
      };

#ifdef FEAT_EICKT
      extern template float DotProduct<Mem::Main>::value_generic(const float * const, const float * const, const Index);
      extern template double DotProduct<Mem::Main>::value_generic(const double * const, const double * const, const Index);
#endif

      template <>
      struct DotProduct<Mem::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);
      };

      template <typename Mem_>
      struct TripleDotProduct;

      template <>
      struct TripleDotProduct<Mem::Main>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
        {
          return value_generic(x, y, z, size);
        }

#ifdef FEAT_HAVE_MKL
        static float value(const float * const x, const float * const y, const float * const z, const Index size)
        {
          return value_mkl(x, y, z, size);
        }

        static double value(const double * const x, const double * const y, const double * const z, const Index size)
        {
          return value_mkl(x, y, z, size);
        }
#endif // FEAT_HAVE_MKL

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
        static __float128 value(const __float128 * const x, const __float128 * const y, const __float128 * const z, const Index size)
        {
          return value_generic(x, y, z, size);
        }
#endif

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);

        static float value_mkl(const float * const x, const float * const y, const float * const z, const Index size);
        static double value_mkl(const double * const x, const double * const y, const double * const z, const Index size);
      };

#ifdef FEAT_EICKT
      extern template float TripleDotProduct<Mem::Main>::value_generic(const float * const, const float * const, const float * const, const Index);
      extern template double TripleDotProduct<Mem::Main>::value_generic(const double * const, const double * const, const double * const, const Index);
#endif

      template <>
      struct TripleDotProduct<Mem::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/dot_product_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
