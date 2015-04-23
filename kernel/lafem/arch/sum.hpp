#pragma once
#ifndef KERNEL_LAFEM_ARCH_SUM_HPP
#define KERNEL_LAFEM_ARCH_SUM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct Sum;

      template <>
      struct Sum<Mem::Main>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
        {
#ifdef FEAST_BACKENDS_MKL
          value_mkl(r, x, y, size);
#else
          value_generic(r, x, y, size);
#endif
        }

#if defined(FEAST_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void value(__float128 * r, const __float128 * const x, const __float128 * const y, const Index size)
        {
          value_generic(r, x, y, size);
        }
#endif

        template <typename DT_>
        static void value_generic(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);

        static void value_mkl(float * r, const float * const x, const float * const y, const Index size);
        static void value_mkl(double * r, const double * const x, const double * const y, const Index size);
      };

      extern template void Sum<Mem::Main>::value_generic(float *, const float * const, const float * const, const Index);
      extern template void Sum<Mem::Main>::value_generic(double *, const double * const, const double * const, const Index);

      template <>
      struct Sum<Mem::CUDA>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/sum_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SUM_HPP
