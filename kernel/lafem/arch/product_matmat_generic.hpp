// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <iostream>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_>
      void ProductMatMat<Mem::Main>::dense_generic(DT_ * r, const DT_ * const x, const DT_ * const y, const Index rows, const Index columns, const Index inner)
      {
        for (Index i(0) ; i < rows ; ++i)
        {
          for (Index j(0) ; j < columns ; ++j)
          {
            DT_ sum(0);
            Index xindex(i * inner);
            Index yindex(j);
            for (Index xcol(0) ; xcol < inner ; ++xcol)
            {
              sum += x[xindex + xcol] * y[yindex + xcol * columns];
            }
            r[i * columns + j] = sum;
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_PRODUCT_MATMAT_GENERIC_HPP
