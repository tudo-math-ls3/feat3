// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_ROW_NORM_HPP
#define KERNEL_LAFEM_ARCH_ROW_NORM_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template <typename Mem_>
      struct RowNorm;

      template <>
      struct RowNorm<Mem::Main>
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

        template <typename DT_, typename IT_>
        static void csr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth)
        {
          bcsr_generic_norm2(row_norms, val, col_ind, row_ptr, rows, BlockHeight, BlockWidth);
        }

        template <typename DT_, typename IT_>
        static void bcsr_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        template <typename DT_ , typename IT_>
        static void ell_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
          const IT_* const cs, const IT_* const cl, const Index C, const Index rows)
        {
          ell_generic_norm2(row_norms, val, col_ind, cs, cl, C, rows);
        }

        template <typename DT_, typename IT_>
        static void ell_generic_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
          const IT_* const cs, const IT_* const cl, const Index C, const Index rows);

        ////////////////////////////////////////////////////////////////////////
        // row norm2sqr
        ////////////////////////////////////////////////////////////////////////
        template <typename DT_, typename IT_>
        static void csr_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows)
        {
          csr_generic_norm2sqr(row_norms, val, col_ind, row_ptr, rows);
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

        template <typename DT_, typename IT_>
        static void bcsr_generic_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_ * const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);

        template <typename DT_ , typename IT_>
        static void ell_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
          const IT_* const cs, const IT_* const cl, const Index C, const Index rows)
        {
          ell_generic_norm2sqr(row_norms, val, col_ind, cs, cl, C, rows);
        }

        template <typename DT_, typename IT_>
        static void ell_generic_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
          const IT_* const cs, const IT_* const cl, const Index C, const Index rows);

        ////////////////////////////////////////////////////////////////////////
        // scaled row norm2sqr
        ////////////////////////////////////////////////////////////////////////
        template <typename DT_, typename IT_>
        static void csr_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows)
        {
          csr_generic_scaled_norm2sqr(row_norms, scal, val, col_ind, row_ptr, rows);
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

        template <typename DT_, typename IT_>
        static void bcsr_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const row_ptr, const Index rows,
        const int BlockHeight, const int BlockWidth);

        template <typename DT_ , typename IT_>
        static void ell_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const cs, const IT_* const cl, const Index C, const Index rows)
        {
          ell_generic_scaled_norm2sqr(row_norms, scal, val, col_ind, cs, cl, C, rows);
        }

        template <typename DT_, typename IT_>
        static void ell_generic_scaled_norm2sqr(DT_* row_norms, const DT_* const scal, const DT_* const val,
        const IT_* const col_ind, const IT_* const cs, const IT_* const cl, const Index C, const Index rows);
      };

#ifdef FEAT_EICKT
      extern template void RowNorm<Mem::Main>::csr_generic_norm2(float*,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm<Mem::Main>::csr_generic_norm2(double*,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm<Mem::Main>::bcsr_generic_norm2(float*,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm<Mem::Main>::bcsr_generic_norm2(double*,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);

      extern template void RowNorm<Mem::Main>::ell_generic_norm2(float*,
      const float* const, const Index* const, const Index* const, const Index* const, Index, const Index);
      extern template void RowNorm<Mem::Main>::ell_generic_norm2(double*,
      const double* const, const Index* const, const Index* const, const Index* const, Index, const Index);

      extern template void RowNorm<Mem::Main>::csr_generic_norm2sqr(float*,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm<Mem::Main>::csr_generic_norm2sqr(double*,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm<Mem::Main>::bcsr_generic_norm2sqr(float*,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm<Mem::Main>::bcsr_generic_norm2sqr(double*,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);

      extern template void RowNorm<Mem::Main>::ell_generic_norm2sqr(float*,
      const float* const, const Index* const, const Index* const, const Index* const, Index, const Index);
      extern template void RowNorm<Mem::Main>::ell_generic_norm2sqr(double*,
      const double* const, const Index* const, const Index* const, const Index* const, Index, const Index);

      extern template void RowNorm<Mem::Main>::csr_generic_scaled_norm2sqr(float*, const float* const,
      const float* const, const Index* const, const Index * const, const Index);
      extern template void RowNorm<Mem::Main>::csr_generic_scaled_norm2sqr(double*, const double* const,
      const double* const, const Index* const, const Index* const, const Index);

      extern template void RowNorm<Mem::Main>::bcsr_generic_scaled_norm2sqr(float*, const float* const,
      const float* const, const Index* const, const Index* const, const Index, const int, const int);
      extern template void RowNorm<Mem::Main>::bcsr_generic_scaled_norm2sqr(double*, const double* const,
      const double* const, const Index* const, const Index* const, const Index, const int, const int);

      extern template void RowNorm<Mem::Main>::ell_generic_scaled_norm2sqr(float*, const float* const,
      const float* const, const Index* const, const Index* const, const Index* const, Index, const Index);
      extern template void RowNorm<Mem::Main>::ell_generic_scaled_norm2sqr(double*, const double* const,
      const double* const, const Index* const, const Index* const, const Index* const, Index, const Index);
#endif


      template <>
      struct RowNorm<Mem::CUDA>
      {
        template <typename DT_, typename IT_>
        static void csr_norm2(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);
        template <typename DT_, typename IT_>
        static void csr_norm2sqr(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);
        template <typename DT_, typename IT_>
        static void csr_scaled_norm2sqr(DT_ * row_norms, const DT_ * const scal, const DT_ * const val, const IT_ * const col_ind, const IT_ * const row_ptr, const Index rows);

        template <typename DT_, typename IT_>
        static void ell_norm2(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind,
          const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows);
        template <typename DT_, typename IT_>
        static void ell_norm2sqr(DT_ * row_norms, const DT_ * const val, const IT_ * const col_ind,
          const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows);
        template <typename DT_, typename IT_>
        static void ell_scaled_norm2sqr(DT_ * row_norms, const DT_ * const scal, const DT_ * const val, const IT_ * const col_ind,
          const IT_ * const cs, const IT_ * const cl, const Index C, const Index rows);

        template <typename DT_, typename IT_>
        static void bcsr_norm2(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);
        template <typename DT_, typename IT_>
        static void bcsr_norm2sqr(DT_* row_norms, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);
        template <typename DT_, typename IT_>
        static void bcsr_scaled_norm2sqr(DT_* row_norms, const DT_ * const scal, const DT_* const val, const IT_* const col_ind,
        const IT_* const row_ptr, const Index rows, const int BlockHeight, const int BlockWidth);
      };

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef  __CUDACC__
#include <kernel/lafem/arch/row_norm_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_ROW_NORM_HPP
