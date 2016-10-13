#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <typeinfo>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct ProductMatMat;

      template <>
      struct ProductMatMat<Mem::Main>
      {
        template <typename DT_>
        static void dense(DT_ * r, const DT_ * const x, const DT_ * const y, const Index rows, const Index columns, const Index inner)
        {
#ifdef FEAT_HAVE_MKL
          dense_mkl(r, x, y, rows, columns, inner);
#else
          dense_generic(r, x, y, rows, columns, inner);
#endif
        }

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void dense(__float128 * r, const __float128 * const x, const __float128 * const y, const Index rows, const Index columns, const Index inner)
        {
          dense_generic(r, x, y, rows, columns, inner);
        }
#endif

        template <typename DT_>
        static void dense_generic(DT_ * r, const DT_ * const x, const DT_ * const y, const Index rows, const Index columns, const Index inner);

        static void dense_mkl(float * r, const float * const x, const float * const y, const Index rows, const Index columns, const Index inner);
        static void dense_mkl(double * r, const double * const x, const double * const y, const Index rows, const Index columns, const Index inner);

      };

#ifdef FEAT_EICKT
      extern template void ProductMatMat<Mem::Main>::dense_generic(float *, const float * const, const float * const, const Index, const Index, const Index);
      extern template void ProductMatMat<Mem::Main>::dense_generic(double *, const double * const, const double * const, const Index, const Index, const Index);
#endif

      template <>
      struct ProductMatMat<Mem::CUDA>
      {
        template <typename DT_>
        static void dense(DT_ * r, const DT_ * const x, const DT_ * const y, const Index rows, const Index columns, const Index inner);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/product_matmat_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
