#pragma once
#ifndef KERNEL_LAFEM_ARCH_NORM_HPP
#define KERNEL_LAFEM_ARCH_NORM_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_, typename VectorT_>
      class Norm2GatewayBase
      {
      public:
        virtual typename VectorT_::DataType value(const VectorT_& x) const = 0;

        virtual ~Norm2GatewayBase()
        {
        }
      };

      template <typename Mem_, typename VectorT_>
      class Norm2SquaredGatewayBase
      {
      public:
        virtual typename VectorT_::DataType value(const VectorT_& x) const = 0;

        virtual ~Norm2SquaredGatewayBase()
        {
        }
      };

      template <typename Mem_>
      struct Norm2;

      template <>
      struct Norm2<Mem::Main>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const Index size)
        {
#ifdef FEAST_BACKENDS_MKL
          return value_mkl(x, size);
#else
          return value_generic(x, size);
#endif
        }

#if defined(FEAST_HAVE_QUADMATH) && !defined(__CUDACC__)
        static __float128 value(const __float128 * const x, const Index size)
        {
          return value_generic(x, size);
        }
#endif

        template <typename DT_>
        static DT_ value_generic(const DT_ * const x, const Index size);

        static float value_mkl(const float * const x, const Index size);
        static double value_mkl(const double * const x, const Index size);
      };

      extern template float Norm2<Mem::Main>::value_generic(const float * const, const Index);
      extern template double Norm2<Mem::Main>::value_generic(const double * const, const Index);

      template <>
      struct Norm2<Mem::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/norm_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_NORM2_HPP
