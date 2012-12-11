#pragma once
#ifndef KERNEL_LAFEM_SUM_HPP
#define KERNEL_LAFEM_SUM_HPP 1

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
    template <typename Algo_>
    struct Sum
    {
    };

    /**
     * \brief Summation calculations.
     *
     * This class calculates summation operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Sum<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x + y\f$
       *
       * \param[out] r The summation result.
       * \param[in] x.The first summand.
       * \param[in] y The second summand.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & x, const DenseVector<Mem::Main, DT_> & y)
      {
        if (x.size() != y.size())
          throw InternalError("Vector size does not match!");
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        const DT_ * yp(y.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        if (rp == xp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] += yp[i];
          }
        }
        else if (rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] += xp[i];
          }
        }
        else if (rp == xp && rp == yp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] += rp[i];
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] + yp[i];
          }
        }
      }

      /**
       * \brief Calculate \f$r \leftarrow x + y\f$
       *
       * \param[out] r The summation result.
       * \param[in] x.The first summand.
       * \param[in] y The second summand.
       */
      template <typename DT_>
      static void value(SparseMatrixCSR<Mem::Main, DT_> & r, const SparseMatrixCSR<Mem::Main, DT_> & x, const SparseMatrixCSR<Mem::Main, DT_> & y)
      {
        if (x.rows() != y.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError("Matrix columns do not match!");
        if (x.rows() != r.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != r.columns())
          throw InternalError("Matrix columns do not match!");

        DenseVector<Mem::Main, DT_> xv(x.used_elements(), x.val());
        DenseVector<Mem::Main, DT_> yv(y.used_elements(), y.val());
        DenseVector<Mem::Main, DT_> rv(r.used_elements(), r.val());

        Sum<Algo::Generic>::value(rv, xv, yv);
      }
    };

    template <>
    struct Sum<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y);
    };

    template <>
    struct Sum<Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DenseVector<Mem::CUDA, DT_> & y);

      template <typename DT_>
      static void value(SparseMatrixCSR<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & x, const SparseMatrixCSR<Mem::CUDA, DT_> & y)
      {
        if (x.rows() != y.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError("Matrix columns do not match!");
        if (x.rows() != r.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != r.columns())
          throw InternalError("Matrix columns do not match!");
        if (x.hash() != r.hash() || y.hash() != r.hash())
          throw InternalError("Matrix layouts do not match!");

        DenseVector<Mem::CUDA, DT_> xv(x.used_elements(), x.val());
        DenseVector<Mem::CUDA, DT_> yv(y.used_elements(), y.val());
        DenseVector<Mem::CUDA, DT_> rv(r.used_elements(), r.val());

        Sum<Algo::CUDA>::value(rv, xv, yv);
      }
    };


  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SUM_HPP
