// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_LUMPING_HPP
#error "Do not include this implementation-only header file directly!"
#endif

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void Lumping::csr_generic(DT_ * lump, const DT_ * const val, const IT_ * const /*col_ind*/,
        const IT_ * const row_ptr, const Index rows)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          DT_ sum(0);
          for (Index col = row_ptr[row]; col < end; col++)
          {
            sum += val[col];
          }
          lump[row] = sum;
        }
      }

      template <typename DT_, typename IT_>
      void Lumping::bcsr_generic(DT_ * lump, const DT_ * const val, const IT_ * const /*col_ind*/,
        const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
      {
        Index block_height = Index(BlockHeight);
        Index block_width = Index(BlockWidth);

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];

          for(Index i(0); i < block_height; ++i)
          {
             lump[block_height*row + i] = DT_(0);
          }

          for (Index col = row_ptr[row]; col < end; col++)
          {
            for(Index i(0); i < block_height; ++i)
            {
              for(Index j(0); j < block_width; ++j)
              {
                lump[block_height*row + i] += val[block_height*block_width*col + i*block_width + j];
              }
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_LUMPING_GENERIC_HPP
