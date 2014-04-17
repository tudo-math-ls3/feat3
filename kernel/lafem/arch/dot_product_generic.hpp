#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
  #error "Do not include this implementation-only header file directly!"
#endif

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template <typename DT_>
DT_ DotProduct<Mem::Main, Algo::Generic>::value(const DT_ * const x, const DT_ * const y, const Index size)
{
  DT_ r(0);

  if(x == y)
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r += x[i] * x[i];
    }
  }
  else
  {
    for (Index i(0) ; i < size ; ++i)
    {
      r += x[i] * y[i];
    }
  }

  return r;
}

#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP
