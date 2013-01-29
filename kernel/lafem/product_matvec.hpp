#pragma once
#ifndef KERNEL_LAFEM_PRODUCT_MATVEC_HPP
#define KERNEL_LAFEM_PRODUCT_MATVEC_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct ProductMatVec
    {
    };

    /**
     * \brief ProductMatVec calculations.
     *
     * This class calculates product operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct ProductMatVec<Algo::Generic>
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

      /**
       * \brief Calculate Matrix-Vector-Product \f$r \leftarrow Ab\f$
       *
       * \param[out] r The product result.
       * \param[in] A The matrix.
       * \param[in] b The vector.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const SparseMatrixELL<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & bv)
      {
        if (bv.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");

        DT_ * result(r.elements());
        const Index * Aj(a.Aj());
        const DT_ * Ax(a.Ax());
        const Index * Arl(a.Arl());
        const DT_ * b(bv.elements());
        const Index stride(a.stride());
        const Index rows(a.rows());

        for (Index row(0) ; row < rows ; ++row)
        {
          const Index * tAj(Aj);
          const DT_ * tAx(Ax);
          DT_ sum(0);
          tAj += row;
          tAx += row;

          const Index max(Arl[row]);
          for(Index n(0); n < max ; n++)
          {
              const DT_ A_ij = *tAx;

              const Index col = *tAj;
              sum += A_ij * b[col];

              tAj += stride;
              tAx += stride;
          }
          result[row] = sum;
        }
      }

      /**
       * \brief Calculate Matrix-Vector-Product \f$r \leftarrow Ab\f$
       *
       * \param[out] r The product result.
       * \param[in] A The matrix.
       * \param[in] b The vector.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const SparseMatrixCOO<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * bp(b.elements());
        const DT_ * val(a.val());
        const Index * row_ptr(a.row());
        const Index * col_ptr(a.column());
        DT_ * rp(r.elements());
        const Index ue(a.used_elements());

        r.clear(DT_(0));
        for (Index i(0) ; i < ue ; ++i)
        {
          rp[row_ptr[i]] += val[i] * bp[col_ptr[i]];
        }
      }
    };

    template <>
    struct ProductMatVec<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const SparseMatrixCSR<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b);
      static void value(DenseVector<Mem::Main, double> & r, const SparseMatrixCSR<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b);

      static void value(DenseVector<Mem::Main, float> & r, const SparseMatrixCOO<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b);
      static void value(DenseVector<Mem::Main, double> & r, const SparseMatrixCOO<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b);
    };

    template <>
    struct ProductMatVec<Algo::CUDA>
    {
      template <typename DT_>
        static void value(DenseVector<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b);

      template <typename DT_>
        static void value(DenseVector<Mem::CUDA, DT_> & r, const SparseMatrixELL<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_PRODUCT_MATVEC_HPP
