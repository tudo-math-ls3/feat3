#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_PRODUCT_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_PRODUCT_HPP 1

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
      struct ComponentProduct;

      template <>
      struct ComponentProduct<Mem::Main, Algo::Generic>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);
      };

      extern template void ComponentProduct<Mem::Main, Algo::Generic>::value(float *, const float * const, const float * const, const Index);
      extern template void ComponentProduct<Mem::Main, Algo::Generic>::value(double *, const double * const, const double * const, const Index);

      template <>
        struct ComponentProduct<Mem::Main, Algo::MKL>
      {
        static void value(float * r, const float * const x, const float * const y, const Index size);
        static void value(double * r, const double * const x, const double * const y, const Index size);
      };

      template <>
      struct ComponentProduct<Mem::CUDA, Algo::CUDA>
      {
        template <typename DT_>
        static void value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#ifndef  __CUDACC__
#include <kernel/lafem/arch/component_product_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_COMPONENT_PRODUCT_HPP
