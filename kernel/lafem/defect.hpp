#pragma once
#ifndef KERNEL_LAFEM_DEFECT_HPP
#define KERNEL_LAFEM_DEFECT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/difference.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct Defect
    {
    };

    /**
     * \brief Defect calculations.
     *
     * This class calculates defect operations.
     *
     * \author Dirk Ribbrock
     */
    template <>
    struct Defect<Algo::Generic>
    {
      /**
       * \brief Calculate \f$r \leftarrow rhs - Ab\f$
       *
       * \param[out] r The defect result.
       * \param[in] rhs The right hand side of the system.
       * \param[in] A The system matrix.
       * \param[in] b The given solution.
         */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & rhs, const SparseMatrixCSR<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");
        if (a.rows() != rhs.size())
          throw InternalError("Vector size does not match!");

        const DT_ * bp(b.elements());
        const DT_ * rhsp(rhs.elements());
        const Index * col_ind(a.col_ind());
        const DT_ * val(a.val());
        const Index * row_ptr(a.row_ptr());
        DT_ * rp(r.elements());
        const Index rows(a.rows());

        for (Index row(0) ; row < rows ; ++row)
        {
          DT_ sum(0);
          const Index end(row_ptr[row + 1]);
          for (Index i(row_ptr[row]) ; i < end ; ++i)
          {
            sum += val[i] * bp[col_ind[i]];
          }
          rp[row] = rhsp[row] - sum;
        }
      }

      /**
       * \brief Calculate \f$r \leftarrow rhs - Ab\f$
       *
       * \param[out] r The defect result.
       * \param[in] rhs The right hand side of the system.
       * \param[in] A The system matrix.
       * \param[in] b The given solution.
         */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & rhs, const SparseMatrixELL<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");
        if (a.rows() != rhs.size())
          throw InternalError("Vector size does not match!");

        DT_ * result(r.elements());
        const DT_ * rhsp(rhs.elements());
        const Index * Aj(a.Aj());
        const DT_ * Ax(a.Ax());
        const Index * Arl(a.Arl());
        const DT_ * bp(b.elements());
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
              sum += A_ij * bp[col];

              tAj += stride;
              tAx += stride;
          }
          result[row] = rhsp[row] - sum;
        }
      }

      /**
       * \brief Calculate \f$r \leftarrow rhs - Ab\f$
       *
       * \param[out] r The defect result.
       * \param[in] rhs The right hand side of the system.
       * \param[in] A The system matrix.
       * \param[in] b The given solution.
         */
      template <typename DT_>
      static void value(DenseVector<Mem::Main, DT_> & r, const DenseVector<Mem::Main, DT_> & rhs, const SparseMatrixCOO<Mem::Main, DT_> & a, const DenseVector<Mem::Main, DT_> & b)
      {
        if (b.size() != a.columns())
          throw InternalError("Vector size does not match!");
        if (a.rows() != r.size())
          throw InternalError("Vector size does not match!");
        if (a.rows() != rhs.size())
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

        Difference<Algo::Generic>::value(r, rhs, r);

      }
    };

    template <>
    struct Defect<Algo::MKL>
    {
      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & rhs, const SparseMatrixCSR<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & rhs, const SparseMatrixCSR<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b);

      static void value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & rhs, const SparseMatrixCOO<Mem::Main, float> & a, const DenseVector<Mem::Main, float> & b);
      static void value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & rhs, const SparseMatrixCOO<Mem::Main, double> & a, const DenseVector<Mem::Main, double> & b);
    };

    template <>
    struct Defect<Algo::CUDA>
    {
      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & rhs, const SparseMatrixCSR<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b);

      template <typename DT_>
      static void value(DenseVector<Mem::CUDA, DT_> & r, const DenseVector<Mem::CUDA, DT_> & rhs, const SparseMatrixELL<Mem::CUDA, DT_> & a, const DenseVector<Mem::CUDA, DT_> & b);
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DEFECT_HPP
