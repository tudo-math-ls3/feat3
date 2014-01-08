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
      template <typename Mem_, typename Algo_>
      struct DotProduct
      {
      };

      template <>
      struct DotProduct<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static DT_ value(const DT_ * const x, const DT_ * const y, const Index size);
      };

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
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
