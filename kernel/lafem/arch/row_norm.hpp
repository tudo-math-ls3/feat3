// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_ROW_NORM_HPP
#define KERNEL_LAFEM_ARCH_ROW_NORM_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct RowNorm
      {
        ////////////////////////////////////////////////////////////////////////
        // row norm2
        ////////////////////////////////////////////////////////////////////////
        template <typename DT_, typename IT_>
        static void csr_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows)
        {
          csr_generic_norm2(row_norms, val, col_ind, row_ptr, rows);
        }

        static void csr_norm2(float* row_norms, const float* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2, csr_generic_norm2, csr_generic_norm2, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2(double* row_norms, const double* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2, csr_generic_norm2, csr_generic_norm2, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2(float* row_norms, const float* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2, csr_generic_norm2, csr_generic_norm2, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2(double* row_norms, const double* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2, csr_generic_norm2, csr_generic_norm2, row_norms, val, col_ind, row_ptr, rows)
        }

        template <typename DT_, typename IT_>
        static void csr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          bcsr_generic_norm2(row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
        }

        static void bcsr_norm2(float* row_norms, const float* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2, bcsr_generic_norm2, bcsr_generic_norm2, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2(double* row_norms, const double* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2, bcsr_generic_norm2, bcsr_generic_norm2, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2(float* row_norms, const float* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2, bcsr_generic_norm2, bcsr_generic_norm2, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2(double* row_norms, const double* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2, bcsr_generic_norm2, bcsr_generic_norm2, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        template <typename DT_, typename IT_>
        static void bcsr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        template <typename DT_, typename IT_>
        static void csr_cuda_norm2(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);
        template <typename DT_, typename IT_>
        static void bcsr_cuda_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        ////////////////////////////////////////////////////////////////////////
        // row norm2sqr
        ////////////////////////////////////////////////////////////////////////
        template <typename DT_, typename IT_>
        static void csr_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows)
        {
          csr_generic_norm2sqr(row_norms, val, col_ind, row_ptr, rows);
        }

        static void csr_norm2sqr(float* row_norms, const float* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2sqr, csr_generic_norm2sqr, csr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2sqr(double* row_norms, const double* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2sqr, csr_generic_norm2sqr, csr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2sqr(float* row_norms, const float* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2sqr, csr_generic_norm2sqr, csr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows)
        }

        static void csr_norm2sqr(double* row_norms, const double* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_norm2sqr, csr_generic_norm2sqr, csr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows)
        }

        template <typename DT_, typename IT_>
        static void csr_generic_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          bcsr_generic_norm2sqr(row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
        }

        static void bcsr_norm2sqr(float* row_norms, const float* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2sqr, bcsr_generic_norm2sqr, bcsr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2sqr(double* row_norms, const double* const val, const unsigned long* const col_ind,
        const unsigned long* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2sqr, bcsr_generic_norm2sqr, bcsr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2sqr(float* row_norms, const float* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2sqr, bcsr_generic_norm2sqr, bcsr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_norm2sqr(double* row_norms, const double* const val, const unsigned int* const col_ind,
        const unsigned int* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_norm2sqr, bcsr_generic_norm2sqr, bcsr_generic_norm2sqr, row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        template <typename DT_, typename IT_>
        static void bcsr_generic_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        template <typename DT_, typename IT_>
        static void csr_cuda_norm2sqr(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);
        template <typename DT_, typename IT_>
        static void bcsr_cuda_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        ////////////////////////////////////////////////////////////////////////
        // scaled row norm2sqr
        ////////////////////////////////////////////////////////////////////////
        template <typename DT_, typename IT_>
        static void csr_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows)
        {
          csr_generic_scaled_norm2sqr(row_norms, scal, val, col_ind, row_ptr, rows);
        }

        static void csr_scaled_norm2sqr(float* row_norms, const float* const scal, const float* const val,
        const unsigned long* const col_ind, const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_scaled_norm2sqr, csr_generic_scaled_norm2sqr, csr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows)
        }

        static void csr_scaled_norm2sqr(double* row_norms, const double* const scal, const double* const val,
        const unsigned long* const col_ind, const unsigned long* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_scaled_norm2sqr, csr_generic_scaled_norm2sqr, csr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows)
        }

        static void csr_scaled_norm2sqr(float* row_norms, const float* const scal, const float* const val,
        const unsigned int* const col_ind, const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_scaled_norm2sqr, csr_generic_scaled_norm2sqr, csr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows)
        }

        static void csr_scaled_norm2sqr(double* row_norms, const double* const scal, const double* const val,
        const unsigned int* const col_ind, const unsigned int* const row_ptr, const Index rows)
        {
          BACKEND_SKELETON_VOID(csr_cuda_scaled_norm2sqr, csr_generic_scaled_norm2sqr, csr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows)
        }

        template <typename DT_, typename IT_>
        static void csr_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth)
        {
          bcsr_generic_scaled_norm2sqr(row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
        }

        static void bcsr_scaled_norm2sqr(float* row_norms, const float* const scal, const float* const val,
        const unsigned long* const col_ind, const unsigned long* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_scaled_norm2sqr(double* row_norms, const double* const scal, const double* const val,
        const unsigned long* const col_ind, const unsigned long* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_scaled_norm2sqr(float* row_norms, const float* const scal, const float* const val,
        const unsigned int* const col_ind, const unsigned int* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        static void bcsr_scaled_norm2sqr(double* row_norms, const double* const scal, const double* const val,
        const unsigned int* const col_ind, const unsigned int* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth)
        {
          BACKEND_SKELETON_VOID(bcsr_cuda_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, bcsr_generic_scaled_norm2sqr, row_norms, scal, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth)
        }

        template <typename DT_, typename IT_>
        static void bcsr_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth);


        template <typename DT_, typename IT_>
        static void csr_cuda_scaled_norm2sqr(DT_ * row_norms, const DT_ * const scal, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);
        template <typename DT_, typename IT_>
        static void bcsr_cuda_scaled_norm2sqr(DT_* row_norms, const DT_ * const scal, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);
      };

#ifdef FEAT_EICKT
      extern template void RowNorm::csr_generic_norm2(float*,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm::csr_generic_norm2(double*,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm::bcsr_generic_norm2(float*,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm::bcsr_generic_norm2(double*,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);

      extern template void RowNorm::csr_generic_norm2sqr(float*,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm::csr_generic_norm2sqr(double*,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm::bcsr_generic_norm2sqr(float*,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm::bcsr_generic_norm2sqr(double*,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);

      extern template void RowNorm::csr_generic_scaled_norm2sqr(float*, const float* const,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm::csr_generic_scaled_norm2sqr(double*, const double* const,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm::bcsr_generic_scaled_norm2sqr(float*, const float* const,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm::bcsr_generic_scaled_norm2sqr(double*, const double* const,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);
#endif


    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/row_norm_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_ROW_NORM_HPP
