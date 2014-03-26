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
      template <typename Mem_, typename Algo_, typename VectorT_>
      class Norm2GatewayBase
      {
        public:
          virtual typename VectorT_::DataType value(const VectorT_& x) const = 0;

          virtual ~Norm2GatewayBase()
          {
          }
      };

      template <typename Mem_, typename Algo_, typename VectorT_>
      class Norm2SquaredGatewayBase
      {
        public:
          virtual typename VectorT_::DataType value(const VectorT_& x) const = 0;

          virtual ~Norm2SquaredGatewayBase()
          {
          }
      };

      template <typename Mem_, typename Algo_>
      struct Norm2;

      template <>
      struct Norm2<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const Index size);
      };

      extern template float Norm2<Mem::Main, Algo::Generic>::value(const float * const, const Index);
      extern template double Norm2<Mem::Main, Algo::Generic>::value(const double * const, const Index);

      template <>
      struct Norm2<Mem::Main, Algo::MKL>
      {
        static float value(const float * const x, const Index size);
        static double value(const double * const x, const Index size);
      };

      template <>
      struct Norm2<Mem::CUDA, Algo::CUDA>
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
