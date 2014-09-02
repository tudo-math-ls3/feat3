#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP 1

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
      class DotGatewayBase
      {
        public:
          virtual typename VectorT_::DataType value(const VectorT_& x, const VectorT_& y) const = 0;

          virtual ~DotGatewayBase()
          {
          }
      };

      // Dot Product
      template <typename Mem_, typename Algo_>
      struct DotProduct;

      template <>
      struct DotProduct<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);
      };

      extern template float DotProduct<Mem::Main, Algo::Generic>::value(const float * const, const float * const, const Index);
      extern template double DotProduct<Mem::Main, Algo::Generic>::value(const double * const, const double * const, const Index);

      template <>
      struct DotProduct<Mem::Main, Algo::MKL>
      {
        static float value(const float * const x, const float * const y, const Index size);

        static double value(const double * const x, const double * const y, const Index size);
      };

      template <>
      struct DotProduct<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

      // Triple dot product
      template <typename Mem_, typename Algo_>
      struct TripleDotProduct;

      template <>
      struct TripleDotProduct<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

      extern template float TripleDotProduct<Mem::Main, Algo::Generic>::value(const float * const, const float * const, const float * const, const Index);
      extern template double TripleDotProduct<Mem::Main, Algo::Generic>::value(const double * const, const double * const, const double * const, const Index);

      template <>
      struct TripleDotProduct<Mem::Main, Algo::MKL>
      {
        static float value(const float * const x, const float * const y, const float * const z, const Index size);
        static double value(const double * const x, const double * const y, const double * const z, const Index size);
      };

      template <>
      struct TripleDotProduct<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

      // Triple dot product, inverted (first vector's components go to the denominator)
      template <typename Mem_, typename Algo_>
      struct TripleDotProductI;

      template <>
      struct TripleDotProductI<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

      extern template float TripleDotProductI<Mem::Main, Algo::Generic>::value(const float * const, const float * const, const float * const, const Index);
      extern template double TripleDotProductI<Mem::Main, Algo::Generic>::value(const double * const, const double * const, const double * const, const Index);

      template <>
      struct TripleDotProductI<Mem::Main, Algo::MKL>
      {
        static float value(const float * const x, const float * const y, const float * const z, const Index size);
        static double value(const double * const x, const double * const y, const double * const z, const Index size);
      };

      template <>
      struct TripleDotProductI<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);

        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/dot_product_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
