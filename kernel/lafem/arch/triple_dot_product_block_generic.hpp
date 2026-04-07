// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_TRIPLE_DOT_PRODUCT_BLOCK_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  template <typename ValueType_>
  ValueType_ TripleDotProductBlock::exec_generic_impl(const ValueType_ * const x, const ValueType_ * const y, const ValueType_ * const z, const Index size)
  {
    ValueType_ r(0);

    if (x == y)
    {

      for(Index i(0); i < size; ++i)
      {
        for(int j(0); j < ValueType_::n; ++j)
        {
          r[j] += x[i][j] * x[i][j] * z[i][j];
        }
      }
    }
    else if (x == z)
    {

      for (Index i(0) ; i < size ; ++i)
      {
        for(int j(0); j < ValueType_::n; ++j)
        {
          r[j] += x[i][j] * x[i][j] * y[i][j];
        }
      }
    }
    else if (y == z)
    {

      for (Index i(0) ; i < size ; ++i)
      {
        for(int j(0); j < ValueType_::n; ++j)
        {
          r[j] += x[i][j] * y[i][j] * y[i][j];
        }
      }
    }
    else
    {

      for(Index i(0); i < size; ++i)
      {
        for(int j(0); j < ValueType_::n; ++j)
        {
          r[j] += x[i][j] * y[i][j] * z[i][j];
        }
      }
    }

    return r;
  }
} // namespace FEAT::LAFEM::Arch
