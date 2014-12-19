#pragma once
#ifndef KERNEL_LAFEM_ARCH_DIFFERENCE_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DIFFERENCE_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DIFFERENCE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {

      template <typename DT_>
      void Difference<Mem::Main, Algo::Generic>::value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
      {
        if (x == r)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] -= y[i];
          }
        }
        else if(y == r)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] = -r[i];
            r[i] += x[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] = x[i] - y[i];
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_DIFFERENCE_GENERIC_HPP
