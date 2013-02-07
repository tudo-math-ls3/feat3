#pragma once
#ifndef KERNEL_LAFEM_TRANSPOSITION_HPP
#define KERNEL_LAFEM_TRANSPOSITION_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>



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
      static SparseMatrixCOO<Mem::Main, DT_> value(const SparseMatrixCOO<Mem::Main, DT_> & a)
      {
        SparseMatrixCOO<Mem::Main, DT_> st(a.columns(), a.rows());

        for (Index i(0) ; i < a.used_elements() ; ++i)
        {
          st(a.column()[i], a.row()[i], a.val()[i]);
        }

        return st;
      }

      template <typename DT_>
      static SparseMatrixCSR<Mem::Main, DT_> value(const SparseMatrixCSR<Mem::Main, DT_> & a)
      {
        SparseMatrixCOO<Mem::Main, DT_> st(a.columns(), a.rows());
        SparseMatrixCOO<Mem::Main, DT_> at(a);

        for (Index i(0) ; i < at.used_elements() ; ++i)
        {
          st(at.column()[i], at.row()[i], at.val()[i]);
        }

        SparseMatrixCSR<Mem::Main, DT_> t(st);
        return t;
      }

      template <typename DT_>
      static SparseMatrixELL<Mem::Main, DT_> value(const SparseMatrixELL<Mem::Main, DT_> & a)
      {
        SparseMatrixCOO<Mem::Main, DT_> st(a.columns(), a.rows());
        SparseMatrixCOO<Mem::Main, DT_> at(a);

        for (Index i(0) ; i < at.used_elements() ; ++i)
        {
          st(at.column()[i], at.row()[i], at.val()[i]);
        }

        SparseMatrixELL<Mem::Main, DT_> t(st);
        return t;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SUM_HPP
