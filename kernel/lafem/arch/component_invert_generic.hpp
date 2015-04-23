#pragma once
#ifndef KERNEL_LAFEM_ARCH_COMPONENT_INVERT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_COMPONENT_INVERT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template<typename DT_>
      void ComponentInvert<Mem::Main>::value_generic(DT_* r, const DT_* const x, const DT_ s, const Index size)
      {
        if (r == x)
        {
          for(Index i(0); i < size; ++i)
          {
            r[i] = s / r[i];
          }
        }
        else
        {
          for(Index i(0); i < size; ++i)
          {
            r[i] = s / x[i];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_COMPONENT_INVERT_HPP
