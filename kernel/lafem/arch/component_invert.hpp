#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template<typename Mem_>
      struct ComponentInvert;

      template<>
      struct ComponentInvert<Mem::Main>
      {
        template<typename DT_>
        static void value(DT_* r, const DT_* const x, const DT_ s, const Index size)
        {
          value_generic(r, x, s, size);
        }

        template<typename DT_>
        static void value_generic(DT_* r, const DT_* const x, const DT_ s, const Index size);
      };

#ifdef FEAT_EICKT
      extern template void ComponentInvert<Mem::Main>::value_generic(float*, const float* const, const float, const Index);
      extern template void ComponentInvert<Mem::Main>::value_generic(double*, const double* const, const double, const Index);
#endif

      template<>
      struct ComponentInvert<Mem::CUDA>
      {
        template<typename DT_>
        static void value(DT_* r, const DT_* const x, const DT_ s, const Index size);
      };
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/component_invert_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP
