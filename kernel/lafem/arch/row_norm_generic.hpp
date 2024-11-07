// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_ROW_NORM_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_ROW_NORM_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_ROW_NORM_HPP
#error "Do not include this implementation-only header file directly!"
#endif

#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename DT_, typename IT_>
      void RowNorm::csr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const /*col_ind*/,
        const IT_* const row_ptr, const Index rows)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          DT_ norm(0);

          // Manually compute norm2 of row
          for (Index col = row_ptr[row]; col < end; col++)
          {
            norm += Math::sqr(val[col]);
          }

          // Take the square root
          row_norms[row] = Math::sqrt(norm);
        }
      }

      template <typename DT_, typename IT_>
      void RowNorm::csr_generic_norm2sqr(DT_* row_norms, const DT_* const val,
      const IT_* const /*col_ind*/, const IT_* const row_ptr, const Index rows)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          DT_ norm(0);

          // Manually compute norm2 of row
          for (Index col = row_ptr[row]; col < end; col++)
          {
            norm += Math::sqr(val[col]);
          }
          // Do not take the square root
          row_norms[row] = norm;
        }
      }

      template <typename DT_, typename IT_>
      void RowNorm::csr_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal,
      const DT_* const val, const IT_* const /*col_ind*/, const IT_* const row_ptr, const Index rows)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];
          DT_ norm(0);

          // Manually compute norm2 of row
          for (Index col = row_ptr[row]; col < end; col++)
          {
            norm += scal[row]*Math::sqr(val[col]);
          }
          // Do not take the square root
          row_norms[row] = norm;
        }
      }

      template <typename DT_, typename IT_>
      void RowNorm::bcsr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const /*col_ind*/,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
      {
        Index block_height = Index(BlockHeight);
        Index block_width = Index(BlockWidth);

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];

          for(Index i = 0; i < block_height; ++i)
          {
             row_norms[block_height*row + i] = DT_(0);
          }

          for (Index col = row_ptr[row]; col < end; col++)
          {
            // Manually compute norm2 of row
            for(Index i = 0; i < block_height; ++i)
            {
              for(Index j = 0; j < block_width; ++j)
              {
                row_norms[block_height*row + i] += Math::sqr(val[block_height*block_width*col + i*block_width + j]);
              }
              // Take the square root
              row_norms[block_height*row + i] = Math::sqrt(row_norms[block_height*row + i]);
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      void RowNorm::bcsr_generic_norm2sqr(DT_* row_norms, const DT_* const val,
      const IT_* const /*col_ind*/, const IT_* const row_ptr, const Index rows,
      const int BlockHeight, const int BlockWidth)
      {
        Index block_height = Index(BlockHeight);
        Index block_width = Index(BlockWidth);

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];

          for(Index i = 0; i < block_height; ++i)
          {
             row_norms[block_height*row + i] = DT_(0);
          }

          for (Index col = row_ptr[row]; col < end; col++)
          {
            // Manually compute norm2 of row
            for(Index i = 0; i < block_height; ++i)
            {
              for(Index j = 0; j < block_width; ++j)
              {
                row_norms[block_height*row + i] += Math::sqr(val[block_height*block_width*col + i*block_width + j]);
              }
              // Do not take the square root
            }
          }
        }
      }

      template <typename DT_, typename IT_>
      void RowNorm::bcsr_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal,
      const DT_* const val, const IT_* const col_ind, const IT_* const row_ptr, const Index rows,
      const int BlockHeight, const int BlockWidth)
      {
        Index block_height = Index(BlockHeight);
        Index block_width = Index(BlockWidth);

        FEAT_PRAGMA_OMP(parallel for)
        for (Index row = 0; row < rows; row++)
        {
          Index end = row_ptr[row + 1];

          for(Index i = 0; i < block_height; ++i)
          {
             row_norms[block_height*row + i] = DT_(0);
          }

          for (Index col = row_ptr[row]; col < end; col++)
          {
            // Manually compute norm2 of row
            for(Index i = 0; i < block_height; ++i)
            {
              for(Index j = 0; j < block_width; ++j)
              {
                row_norms[block_height*row + i] += scal[block_width*col_ind[col] + j]
                  *Math::sqr(val[block_height*block_width*col + i*block_width + j]);
              }
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_ARCH_ROW_NORM_GENERIC_HPP
