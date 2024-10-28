// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_DIAGONAL_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_DIAGONAL_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_DIAGONAL_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename IT_>
      void Diagonal::csr_generic(IT_ * diag, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          const Index end = row_ptr[row + 1];
          diag[row] = row_ptr[rows];
          for (Index col = row_ptr[row]; col < end; col++)
          {
            if (row == col_ind[col])
            {
              diag[row] = IT_(col);
              break;
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_DIAGONAL_GENERIC_HPP
