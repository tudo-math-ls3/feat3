// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DOT_PRODUCT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      DT_ DotProduct<Mem::Main>::value_generic(const DT_ * const x, const DT_ * const y, const Index size)
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

      template <typename DT_>
      DT_ TripleDotProduct<Mem::Main>::value_generic(const DT_ * const x, const DT_ * const y, const DT_ * const z, const Index size)
      {
        DT_ r(0);

        if (x == y)
        {
          for (Index i(0) ; i < size ; ++i)
            r += x[i] * x[i] * z[i];
        }
        if (x == z)
        {
          for (Index i(0) ; i < size ; ++i)
            r += x[i] * x[i] * y[i];
        }
        if (y == z)
        {
          for (Index i(0) ; i < size ; ++i)
            r += x[i] * y[i] * y[i];
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
            r += x[i] * y[i] * z[i];
        }

        return r;
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_DOT_PRODUCT_GENERIC_HPP
