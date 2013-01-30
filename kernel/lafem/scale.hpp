#pragma once
#ifndef KERNEL_LAFEM_SCALE_HPP
#define KERNEL_LAFEM_SCALE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Scale
    {
    };

    /**
     * \brief Scaling calculations.
     *
     * This class calculates scaling operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Scale<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow x \cdot s\f$
       *
       * \param[out] r The scaled vector.
       * \param[in] x The vector to be scaled.
       * \param[in] s A scalar to scale x with.
       */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & x, const DT_ s)
      {
        if (x.size() != r.size())
          throw InternalError("Vector size does not match!");

        const DT_ * xp(x.elements());
        DT_ * rp(r.elements());
        const Index size(r.size());

        if (xp == rp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= s;
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] * s;
          }
        }
      }

      template <typename DT_>
      static void value(SparseMatrixCOO<Mem::Main, DT_> & r, const SparseMatrixCOO<Mem::Main, DT_> & x, const DT_ s)
      {
        if(x.rows() != r.rows())
          throw InternalError("Matrix Rows doe not match!");
        if(x.columns() != r.columns())
          throw InternalError("Matrix Columns doe not match!");

        const DT_ * xp(x.val());
        DT_ * rp(r.val());
        const Index size(r.used_elements());

        if (xp == rp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= s;
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] * s;
          }
        }
      }

      template <typename DT_>
      static void value(SparseMatrixCSR<Mem::Main, DT_> & r, const SparseMatrixCSR<Mem::Main, DT_> & x, const DT_ s)
      {
        if(x.rows() != r.rows())
          throw InternalError("Matrix Rows doe not match!");
        if(x.columns() != r.columns())
          throw InternalError("Matrix Columns doe not match!");

        const DT_ * xp(x.val());
        DT_ * rp(r.val());
        const Index size(r.used_elements());

        if (xp == rp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= s;
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] * s;
          }
        }
      }

      template <typename DT_>
      static void value(SparseMatrixELL<Mem::Main, DT_> & r, const SparseMatrixELL<Mem::Main, DT_> & x, const DT_ s)
      {
        if(x.rows() != r.rows())
          throw InternalError("Matrix Rows doe not match!");
        if(x.columns() != r.columns())
          throw InternalError("Matrix Columns doe not match!");

        const DT_ * xp(x.Ax());
        DT_ * rp(r.Ax());
        const Index size(r.stride() * r.num_cols_per_row());

        if (xp == rp)
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] *= s;
          }
        }
        else
        {
          for (Index i(0) ; i < size ; ++i)
          {
            rp[i] = xp[i] * s;
          }
        }
      }
    };

    template <>
    struct Scale<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const float s);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const double s);

      static void value(SparseMatrixCOO<Mem::Main, float> & r, const SparseMatrixCOO<Mem::Main, float> & x, const float s);
      static void value(SparseMatrixCOO<Mem::Main, double> & r, const SparseMatrixCOO<Mem::Main, double> & x, const double s);

      static void value(SparseMatrixCSR<Mem::Main, float> & r, const SparseMatrixCSR<Mem::Main, float> & x, const float s);
      static void value(SparseMatrixCSR<Mem::Main, double> & r, const SparseMatrixCSR<Mem::Main, double> & x, const double s);

      static void value(SparseMatrixELL<Mem::Main, float> & r, const SparseMatrixELL<Mem::Main, float> & x, const float s);
      static void value(SparseMatrixELL<Mem::Main, double> & r, const SparseMatrixELL<Mem::Main, double> & x, const double s);
    };

    template <>
    struct Scale<Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & x, const DT_ s);

      template <typename DT_>
      static void value(SparseMatrixCOO<Mem::CUDA, DT_> & r, const SparseMatrixCOO<Mem::CUDA, DT_> & x, const DT_ s);

      template <typename DT_>
      static void value(SparseMatrixCSR<Mem::CUDA, DT_> & r, const SparseMatrixCSR<Mem::CUDA, DT_> & x, const DT_ s);

      template <typename DT_>
      static void value(SparseMatrixELL<Mem::CUDA, DT_> & r, const SparseMatrixELL<Mem::CUDA, DT_> & x, const DT_ s);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SCALE_HPP
