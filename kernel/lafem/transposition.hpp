#pragma once
#ifndef KERNEL_LAFEM_TRANSPOSITION_HPP
#define KERNEL_LAFEM_TRANSPOSITION_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Transposition
    {
    };

    /**
     * \brief Summation calculations.
     *
     * This class transposes a given matrix ex-situ.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Transposition<Algo::Generic>
    {
      /**
       * \brief Calculate \f$A^T\f$
       *
       * \param[in] a The input matrix.
       * \returns The transpose of a.
       */
      template <typename DT_>
      static SparseMatrixCSR<Mem::Main, DT_> value(const SparseMatrixCSR<Mem::Main, DT_> & a)
      {
        /*DenseVector<Mem::Main, Index> col_ind(a.used_elements());
        DenseVector<Mem::Main, DT_> val(a.used_elements());
        DenseVector<Mem::Main, Index> row_ptr(a.columns() + 1);
        DenseVector<Mem::Main, Index> row_ptr_end(a.columns());*/

        SparseMatrixCOO<Mem::Main, DT_> st(a.columns(), a.rows());

        const Index * col_ind(a.col_ind());
        const DT_ * val(a.val());
        const Index * row_ptr(a.row_ptr());
        const Index * row_ptr_end(a.row_ptr_end());
        const Index rows(a.rows());

        for (Index row(0) ; row < rows ; ++row)
        {
          const Index end(row_ptr_end[row]);
          for (Index i(row_ptr[row]) ; i < end ; ++i)
          {
            st(col_ind[i], row, val[i]);
          }
        }

        SparseMatrixCSR<Mem::Main, DT_> t(st);
        //SparseMatrixCSR<Mem::Main, DT_> t(a.columns(), a.rows(), col_ind, val, row_ptr, row_ptr_end);
        return t;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SUM_HPP
