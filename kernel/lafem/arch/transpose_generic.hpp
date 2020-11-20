// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_TRANSPOSE_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_TRANSPOSE_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_TRANSPOSE_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      /// \todo cache blocking https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
      template <typename DT_>
      void Transpose<Mem::Main>::value_generic(DT_ * r, const DT_ * const x, const Index rows_x, const Index columns_x)
      {
        for (Index i(0) ; i < rows_x ; ++i)
        {
          for (Index j(0) ; j < columns_x ; ++j)
          {
            r[j * rows_x + i] = x[i * columns_x + j];
          }
        }
      }

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_TRANSPOSE_GENERIC_HPP
