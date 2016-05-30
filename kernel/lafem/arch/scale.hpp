#pragma once
#ifndef KERNEL_LAFEM_ARCH_SCALE_HPP
#define KERNEL_LAFEM_ARCH_SCALE_HPP 1

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
      struct Scale;

      template <>
      struct Scale<Mem::Main>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ s, const Index size)
        {
#ifdef FEAT_BACKENDS_MKL
          value_mkl(r, x, s, size);
#else
          value_generic(r, x, s, size);
#endif
        }

#if defined(FEAT_HAVE_QUADMATH) && !defined(__CUDACC__)
        static void value(__float128 * r, const __float128 * const x, const __float128 s, const Index size)
        {
          value_generic(r, x, s, size);
        }
#endif

        template <typename DT_>
        static void value_generic(DT_ * r, const DT_ * const x, const DT_ s, const Index size);

        static void value_mkl(float * r, const float * const x, const float, const Index size);
        static void value_mkl(double * r, const double * const x, const double, const Index size);
      };

      extern template void Scale<Mem::Main>::value_generic(float *, const float * const, const float, const Index);
      extern template void Scale<Mem::Main>::value_generic(double *, const double * const, const double, const Index);

      template <>
      struct Scale<Mem::CUDA>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/scale_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_SCALE_HPP
