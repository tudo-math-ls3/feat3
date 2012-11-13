#pragma once
#ifndef KERNEL_LAFEM_PRODUCT_HPP
#define KERNEL_LAFEM_PRODUCT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Arch_, typename BType_>
    struct Product
    {
    };

    /**
     * \brief Product calculations.
     *
     * This class calculates product operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Product<Mem::Main, Algo::Generic>
    {
      /**
       * \brief Calculate Matrix-Vector-Product \f$r \leftarrow Ab\f$
       *
       * \param[out] r The product result.
       * \param[in] A The matrix.
       * \param[in] b The vector.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const SparseMatrixCSR<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * bp(b.elements());
        const Index * col_ind(a.col_ind());
        const DT_ * val(a.val());
        const Index * row_ptr(a.row_ptr());
        const Index * row_ptr_end(a.row_ptr_end());
        DT_ * rp(r.elements());
        const Index rows(a.rows());

        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          const Index end(row_ptr_end[row]);
          for (Index i(row_ptr[row]) ; i < end ; ++i)
          {
            sum += val[i] * bp[col_ind[i]];
          }
          rp[row] = sum;
        }
      }
    };

    template <>
    struct Product<Mem::CUDA, Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_PRODUCT_HPP
