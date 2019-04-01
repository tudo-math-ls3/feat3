// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
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
      template <typename DT_, typename IT_>
      void Diagonal<Mem::Main>::csr_generic(DT_ * diag, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
      {
        for (Index row(0); row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          diag[row] = DT_(0);
          for (Index col = row_ptr[row]; col < end; col++)
          {
            if (row == col_ind[col])
            {
              diag[row] = val[col];
              break;
            }
          }
        }
      }

        template <typename DT_, typename IT_, int BlockHeight_, int BlockWidth_>
      void Diagonal<Mem::Main>::csrb_generic(DT_ * diag, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows)
      {
        for (Index row(0); row < rows; row++)
        {
          for (Index i(0) ; i < BlockHeight_ ; ++i)
          {
            diag[row * BlockHeight_ + i] = DT_(0);
          }
          const Index end = row_ptr[row + 1];

          for (Index col = row_ptr[row]; col < end; col++)
          {
            if (row == col_ind[col])
            {
              for (Index i(0) ; i < BlockHeight_ ; ++i)
              {
                diag[row * BlockHeight_ + i] = val[(col * BlockHeight_ * BlockWidth_) + i + i * BlockWidth_];
              }
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      void Diagonal<Mem::Main>::ell_generic(DT_ * diag, const DT_ * const val, const IT_ * const col_ind,
        const IT_ * const cs, const IT_ * const /*cl*/, const Index C, const Index rows)
      {
        for (Index row(0) ; row < rows ; ++row)
        {
          const Index chunk(row / C);
          const Index local_row(row % C);
          const Index chunk_end(cs[chunk+1]);

          diag[row] = DT_(0);
          for (Index pcol(cs[chunk] + local_row) ; pcol < chunk_end ; pcol+=C)
          {
            if (col_ind[pcol] == row)
            {
              diag[row] = val[pcol];
              break;
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_DIAGONAL_GENERIC_HPP
