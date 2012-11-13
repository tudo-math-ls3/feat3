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
    template <typename Arch_, typename BType_>
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
    struct Sum<Archs::CPU, Archs::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x + y\f$
       *
       * \param[out] r The summation result.
       * \param[in] x.The first summand.
       * \param[in] y The second summand.
       */
      template <typename DT_>
      static void value(DenseVector<Archs::CPU, DT_> & r, const DenseVector<Archs::CPU, DT_> & x, const DenseVector<Archs::CPU, DT_> & y)
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
      static void value(SparseMatrixCSR<Archs::CPU, DT_> & r, const SparseMatrixCSR<Archs::CPU, DT_> & x, const SparseMatrixCSR<Archs::CPU, DT_> & y)
      {
        if (x.rows() != y.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError("Matrix columns do not match!");
        if (x.rows() != r.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != r.columns())
          throw InternalError("Matrix columns do not match!");

        DenseVector<Archs::CPU, DT_> xv(x.used_elements(), x.val());
        DenseVector<Archs::CPU, DT_> yv(y.used_elements(), y.val());
        DenseVector<Archs::CPU, DT_> rv(r.used_elements(), r.val());

        Sum<Archs::CPU, Archs::Generic>::value(rv, xv, yv);
      }
    };

    template <>
    struct Sum<Archs::GPU, Archs::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Archs::GPU, DT_> & r, const DenseVector<Archs::GPU, DT_> & x, const DenseVector<Archs::GPU, DT_> & y);

      template <typename DT_>
      static void value(SparseMatrixCSR<Archs::GPU, DT_> & r, const SparseMatrixCSR<Archs::GPU, DT_> & x, const SparseMatrixCSR<Archs::GPU, DT_> & y)
      {
        if (x.rows() != y.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != y.columns())
          throw InternalError("Matrix columns do not match!");
        if (x.rows() != r.rows())
          throw InternalError("Matrix rows do not match!");
        if (x.columns() != r.columns())
          throw InternalError("Matrix columns do not match!");

        DenseVector<Archs::GPU, DT_> xv(x.used_elements(), x.val());
        DenseVector<Archs::GPU, DT_> yv(y.used_elements(), y.val());
        DenseVector<Archs::GPU, DT_> rv(r.used_elements(), r.val());

        Sum<Archs::GPU, Archs::CUDA>::value(rv, xv, yv);
      }
    };


  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SUM_HPP
