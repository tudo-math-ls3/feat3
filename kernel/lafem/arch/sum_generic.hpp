#pragma once
#ifndef KERNEL_LAFEM_ARCH_SUM_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_SUM_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_SUM_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

namespace FEAST
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      void Sum<Mem::Main, Algo::Generic>::value(DT_ * r, const DT_ * const x, const DT_ * const y, const Index size)
      {
        if (r == x)
        {
          for (Index i(0) ; i < size ; ++i)
            r[i] += y[i];
        }
        else if (r == y)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] += x[i];
          }
        }
        else if (r == x && r == y)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] += r[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            r[i] = x[i] + y[i];
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_ARCH_SUM_GENERIC_HPP
