// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef KERNEL_LAFEM_ARCH_TRANSPOSE_DENSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT::LAFEM::Arch
{
  /// \todo cache blocking https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
  template <typename DT_>
  void TransposeDense::exec_generic_impl(DT_ * r, const DT_ * const x, const Index rows_x, const Index cols_x)
  {
    if (r == x)
    {
      /// \todo use inplace transform, i.e. lower triangular swap algorithm
      DT_* t= new DT_[rows_x * cols_x];
      Memory::memcopy_main(t, x, rows_x * cols_x * sizeof(DT_));
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < rows_x ; ++i)
      {
        for (Index j = 0 ; j < cols_x ; ++j)
        {
          r[j * rows_x + i] = t[i * cols_x + j];
        }
      }
      delete[] t;
    }
    else
    {
      FEAT_PRAGMA_OMP(parallel for)
      for (Index i = 0 ; i < rows_x ; ++i)
      {
        for (Index j = 0 ; j < cols_x ; ++j)
        {
          r[j * rows_x + i] = x[i * cols_x + j];
        }
      }
    }
  }
} // namespace FEAT::LAFEM::Arch
