#pragma once
#ifndef KERNEL_LAFEM_BLOCK_TO_SCALAR_HPP
#define KERNEL_LAFEM_BLOCK_TO_SCALAR_HPP 1

// includes, FEAST
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_layout.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    template<typename Algo_>
    struct MatBlockToScalar;

    /**
     * \brief Convert any matrix to scalar matrix
     *
     * This class converts any matrix to a scalar matrix
     *
     * \author Christoph Lohmann
     */
    template<>
    struct MatBlockToScalar<Algo::Generic>
    {
    private:
      typedef Algo::Generic Algo_;

      /**
       * SparseMatrixCSR
       */
      template<typename DT_>
      static Index _get_length_of_line(const SparseMatrixCSR<Mem::Main, DT_> & matrix, const Index row)
      {
        const Index * prow_ptr(matrix.row_ptr());
        return prow_ptr[row + 1] - prow_ptr[row];
      }

      template<typename DT_>
      static void _set_line(const SparseMatrixCSR<Mem::Main, DT_> & matrix,
                            const Index row, DT_ * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index * prow_ptr(matrix.row_ptr());
        const Index * pcol_ind(matrix.col_ind());
        const DT_ * pval(matrix.val());

        const Index start(prow_ptr[row]);
        for (Index i(0); i < prow_ptr[row + 1] - prow_ptr[row]; ++i)
        {
          pval_set[i * stride] = pval[start + i];
          pcol_set[i * stride] = pcol_ind[start + i] + col_start;
        }
      }

      /**
       * SparseMatrixCOO
       */
      template<typename DT_>
      static Index _get_length_of_line(const SparseMatrixCOO<Mem::Main, DT_> & matrix, const Index row)
      {
        const Index * prow(matrix.row_indices());
        const Index used_elements(matrix.used_elements());

        Index i(0);
        while (prow[i] < row)
        {
          ++i;
        }

        Index j(i);
        while (j < used_elements && prow[j] == row)
        {
          ++j;
        }
        return (j - i);
      }

      template<typename DT_>
      static void _set_line(const SparseMatrixCOO<Mem::Main, DT_> & matrix,
                            const Index row, DT_ * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index * prow(matrix.row_indices());
        const Index * pcol(matrix.column_indices());
        const DT_ * pval(matrix.val());

        const Index used_elements(matrix.used_elements());

        Index start(0);
        while (prow[start] < row)
        {
          ++start;
        }

        for (Index i(0); start + i < used_elements; ++i)
        {
          if (prow[start + i] != row)
          {
            return;
          }
          pval_set[i * stride] = pval[start + i];
          pcol_set[i * stride] = pcol[start + i] + col_start;
        }
      }

      /**
       * SparseMatrixELL
       */
      template<typename DT_>
      static Index _get_length_of_line(const SparseMatrixELL<Mem::Main, DT_> & matrix, const Index row)
      {
        return matrix.Arl()[row];
      }

      template<typename DT_>
      static void _set_line(const SparseMatrixELL<Mem::Main, DT_> & matrix,
                            const Index row, DT_ * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index * parl(matrix.Arl());
        const Index * paj(matrix.Aj() + row);
        const DT_ * pax(matrix.Ax() + row);
        const Index astride(matrix.stride());

        const Index length(parl[row]);

        for (Index i(0); i < length; ++i)
        {
          pval_set[i * stride] = pax[i * astride];
          pcol_set[i * stride] = paj[i * astride] + col_start;
        }
      }

      /**
       * SaddlePointMatrix
       */
      template<
        typename MatrixA_,
        typename MatrixB_ = MatrixA_,
        typename MatrixD_ = MatrixB_>
      static Index _get_length_of_line(const SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> & matrix, const Index row)
      {
        const Index arows(matrix.block_a().rows());

        if (row < arows)
        {
          return _get_length_of_line(matrix.block_a(), row) + _get_length_of_line(matrix.block_b(), row);
        }
        else
        {
          return _get_length_of_line(matrix.block_d(), row - arows);
        }
      }

      template<
        typename MatrixA_,
        typename MatrixB_ = MatrixA_,
        typename MatrixD_ = MatrixB_>
      static void _set_line(const SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> & matrix,
                            const Index row, typename MatrixA_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index arows(matrix.block_a().rows());

        if (row < arows)
        {
          const Index length_of_a(_get_length_of_line(matrix.block_a(), row));

          _set_line(matrix.block_a(), row, pval_set, pcol_set, col_start, stride);
          _set_line(matrix.block_b(), row, pval_set + stride * length_of_a, pcol_set + stride * length_of_a, col_start + matrix.block_a().columns(), stride);
        }
        else
        {
          _set_line(matrix.block_d(), row - arows, pval_set, pcol_set, col_start, stride);
        }
      }

      /**
       * DenseMatrix
       */
      template<typename Mem_, typename DT_>
      static Index _get_length_of_line(const DenseMatrix<Mem_, DT_> & matrix, const Index /*row*/)
      {
        return matrix.columns();
      }

      template<typename Mem_, typename DT_>
      static void _set_line(const DenseMatrix<Mem_, DT_> & matrix,
                            const Index row, DT_ * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index acolumns(matrix.columns());
        const Index * pval(matrix.elements() + row * acolumns);

        for (Index i(0); i < acolumns; ++i)
        {
          pval_set[i * stride] = pval[i];
          pcol_set[i * stride] = i + col_start;
        }
      }

      /**
       * PowerColMatrix
       */
      template<typename SubType_, Index blocks_>
      static Index _get_length_of_line(const PowerColMatrix<SubType_, blocks_> & matrix, const Index row)
      {
        const Index brows(matrix.base().rows());

        if (row < brows)
        {
          return _get_length_of_line(matrix.base(), row);
        }
        else
        {
          return _get_length_of_line(matrix.last(), row - brows);
        }
      }

      template<typename SubType_>
      static Index _get_length_of_line(const PowerColMatrix<SubType_, Index(1)> & matrix, const Index row)
      {
        return _get_length_of_line(matrix.last(), row);
      }

      template<typename SubType_, Index blocks_>
      static void _set_line(const PowerColMatrix<SubType_, blocks_> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index brows(matrix.base().rows());

        if (row < brows)
        {
          _set_line(matrix.base(), row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          _set_line(matrix.last(), row - brows, pval_set, pcol_set, col_start, stride);
        }
      }

      template<typename SubType_>
      static void _set_line(const PowerColMatrix<SubType_, Index(1)> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        _set_line(matrix.last(), row, pval_set, pcol_set, col_start, stride);
      }

      /**
       * PowerRowMatrix
       */
      template<typename SubType_, Index blocks_>
      static Index _get_length_of_line(const PowerRowMatrix<SubType_, blocks_> & matrix, const Index row)
      {
        return _get_length_of_line(matrix.base(), row) + _get_length_of_line(matrix.last(), row);
      }

      template<typename SubType_>
      static Index _get_length_of_line(const PowerRowMatrix<SubType_, Index(1)> & matrix, const Index row)
      {
        return _get_length_of_line(matrix.last(), row);
      }

      template<typename SubType_, Index blocks_>
      static void _set_line(const PowerRowMatrix<SubType_, blocks_> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index length_of_base(_get_length_of_line(matrix.base(), row));

        _set_line(matrix.base(), row, pval_set, pcol_set, col_start, stride);
        _set_line(matrix.last(), row, pval_set + stride * length_of_base, pcol_set + stride * length_of_base, col_start + matrix.base().columns(), stride);
      }

      template<typename SubType_>
      static void _set_line(const PowerRowMatrix<SubType_, Index(1)> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        _set_line(matrix.last(), row, pval_set, pcol_set, col_start, stride);
      }

      /**
       * PowerDiagMatrix
       */
      template<typename SubType_, Index blocks_>
      static Index _get_length_of_line(const PowerDiagMatrix<SubType_, blocks_> & matrix, const Index row)
      {
        const Index brows(matrix.base().rows());

        if (row < brows)
        {
          return _get_length_of_line(matrix.base(), row);
        }
        else
        {
          return _get_length_of_line(matrix.base(), row - brows);
        }
      }

      template<typename SubType_>
      static Index _get_length_of_line(const PowerDiagMatrix<SubType_, Index(1)> & matrix, const Index row)
      {
        return _get_length_of_line(matrix.last(), row);
      }

      template<typename SubType_, Index blocks_>
      static void _set_line(const PowerDiagMatrix<SubType_, blocks_> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        const Index brows(matrix.base().rows());
        const Index bcolumns(matrix.base().columns());

        if (row < brows)
        {
          _set_line(matrix.base(), row, pval_set, pcol_set, col_start, stride);
        }
        else
        {
          _set_line(matrix.last(), row - brows, pval_set, pcol_set, col_start + bcolumns, stride);
        }
      }

      template<typename SubType_>
      static void _set_line(const PowerDiagMatrix<SubType_, Index(1)> & matrix,
                            const Index row, typename SubType_::DataType * pval_set, Index * pcol_set, const Index col_start, const Index stride = 1)
      {
        _set_line(matrix.last(), row, pval_set, pcol_set, col_start, stride);
      }

    public:
      /**
       * \brief Converts any matrix to SparseMatrixCSR-format
       *
       * \param[in] a The input matrix.
       * \returns The matrix a in SparseMatrixCSR-format.
       */
      template<typename MT_>
      static SparseMatrixCSR<typename MT_::MemType, typename MT_::DataType> value_csr(const MT_ & a)
      {
        typedef typename MT_::DataType DT_;
        typedef typename MT_::MemType Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index aused_elements(a.used_elements());

        DenseVector<Mem_, DT_> val(aused_elements);
        DenseVector<Mem_, Index> col_ind(aused_elements);
        DenseVector<Mem_, Index> row_ptr(arows + 1);

        DT_ * pval(val.elements());
        Index * pcol_ind(col_ind.elements());
        Index * prow_ptr(row_ptr.elements());

        for (Index i(0); i < arows; ++i)
        {
          prow_ptr[i + 1] = _get_length_of_line(a, i);
        }

        prow_ptr[0] = 0;

        for (Index i(1); i < arows + 1; ++i)
        {
          prow_ptr[i] += prow_ptr[i - 1];
        }

        for (Index i(0); i < arows; ++i)
        {
          _set_line(a, i, pval + prow_ptr[i], pcol_ind + prow_ptr[i], 0);
        }

        SparseMatrixCSR<Mem_, DT_> ta(arows, acolumns, col_ind, val, row_ptr);

        return ta;
      }

      /**
       * \brief Converts any matrix to SparseMatrixCOO-format
       *
       * \param[in] a The input matrix.
       * \returns The matrix a in SparseMatrixCOO-format.
       */
      template<typename MT_>
      static SparseMatrixCOO<typename MT_::MemType, typename MT_::DataType> value_coo(const MT_ & a)
      {
        typedef typename MT_::DataType DT_;
        typedef typename MT_::MemType Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index aused_elements(a.used_elements());

        DenseVector<Mem_, DT_> val(aused_elements);
        DenseVector<Mem_, Index> col_ind(aused_elements);
        DenseVector<Mem_, Index> row_ind(aused_elements);

        DT_ * pval(val.elements());
        Index * pcol_ind(col_ind.elements());
        Index * prow_ind(row_ind.elements());

        DenseVector<Mem_, Index> row_ptr(arows + 1);
        Index * prow_ptr(row_ptr.elements());

        for (Index i(0); i < arows; ++i)
        {
          prow_ptr[i + 1] = _get_length_of_line(a, i);
        }

        prow_ptr[0] = 0;

        for (Index i(1); i < arows + 1; ++i)
        {
          prow_ptr[i] += prow_ptr[i - 1];
        }

        for (Index i(0); i < arows; ++i)
        {
          _set_line(a, i, pval + prow_ptr[i], pcol_ind + prow_ptr[i], 0);
        }

        for (Index i(0); i < arows; ++i)
        {
          for (Index k(prow_ptr[i]); k < prow_ptr[i+1]; ++k)
          {
            prow_ind[k] = i;
          }
        }

        SparseMatrixCOO<Mem_, DT_> ta(arows, acolumns, row_ind, col_ind, val);

        return ta;
      }

      /**
       * \brief Converts any matrix to SparseMatrixELL-format
       *
       * \param[in] a The input matrix.
       * \returns The matrix a in SparseMatrixELL-format.
       */
      template<typename MT_>
      static SparseMatrixELL<typename MT_::MemType, typename MT_::DataType> value_ell(const MT_ & a)
      {
        typedef typename MT_::DataType DT_;
        typedef typename MT_::MemType Mem_;

        const Index arows(a.rows());
        const Index acolumns(a.columns());
        const Index aused_elements(a.used_elements());

        Index alignment(32);
        const Index stride(alignment * ((arows + alignment - 1)/ alignment));

        DenseVector<Mem_, Index> arl(arows);
        Index * parl(arl.elements());

        Index num_cols_per_row(0);
        for (Index i(0); i < arows; ++i)
        {
          parl[i] = _get_length_of_line(a, i);
          if (num_cols_per_row < parl[i])
          {
            num_cols_per_row = parl[i];
          }
        }

        DenseVector<Mem_, DT_> ax(stride * num_cols_per_row);
        DenseVector<Mem_, Index> aj(stride * num_cols_per_row);

        DT_ * pax(ax.elements());
        Index * paj(aj.elements());

        for (Index i(0); i < arows; ++i)
        {
          _set_line(a, i, pax + i, paj + i, 0, stride);
        }

        SparseMatrixELL<Mem_, DT_> ta(arows, acolumns, stride, num_cols_per_row, aused_elements, ax, aj, arl);
        return ta;
      }



      template<SparseLayoutId LT_, typename MT_>
      static SparseMatrixCOO<typename MT_::MemType, typename MT_::DataType> value(const MT_ & a, typename std::enable_if<LT_ == SparseLayoutId::lt_coo>::type* = 0)
      {
        return value_coo(a);
      }

      template<SparseLayoutId LT_, typename MT_>
      static SparseMatrixCSR<typename MT_::MemType, typename MT_::DataType> value(const MT_ & a, typename std::enable_if<LT_ == SparseLayoutId::lt_csr>::type* = 0)
      {
        return value_csr(a);
      }

      template<SparseLayoutId LT_, typename MT_>
      static SparseMatrixELL<typename MT_::MemType, typename MT_::DataType> value(const MT_ & a, typename std::enable_if<LT_ == SparseLayoutId::lt_ell>::type* = 0)
      {
        return value_ell(a);
      }
    };

    template<typename Algo_>
    struct VecBlockToScalar;

    /**
     * \brief Convert any vector to scalar vector
     *
     * This class converts any vector to a scalar vector
     *
     * \author Christoph Lohmann
     */
    template<>
    struct VecBlockToScalar<Algo::Generic>
    {
    private:
      typedef Algo::Generic Algo_;

      /**
       * DenseVector
       */
      template<typename DT_>
      static void _set_vec(const DenseVector<Mem::Main, DT_> & vec, DT_ * pval_set)
      {
        const DT_ * pvec(vec.elements());
        const Index n(vec.size());

        for (Index i(0); i < n; ++i)
        {
          pval_set[i] = pvec[i];
        }
      }

      /**
       * PowerVector
       */
      template<typename VT_, Index count_>
      static void _set_vec(const PowerVector<VT_, count_> & vec, typename VT_::DataType * pval_set)
      {
        _set_vec(vec.base(), pval_set);
        _set_vec(vec.last(), pval_set + vec.base().size());
      }

      template<typename VT_>
      static void _set_vec(const PowerVector<VT_, Index(1)> & vec, typename VT_::DataType * pval_set)
      {
        _set_vec(vec.last(), pval_set);
      }

      /**
       * TupleVector
       */
      template<typename First_, typename Second_, typename... Rest_>
      static void _set_vec(const TupleVector<First_, Second_, Rest_...> & vec, typename TupleVector<First_, Second_, Rest_...>::DataType * pval_set)
      {
        _set_vec(vec.first(), pval_set);
        _set_vec(vec.rest(), pval_set + vec.first().size());
      }

#ifdef FEAST_COMPILER_MICROSOFT
      // compiler hack...
      template<typename... First_>
      static void _set_vec(const TupleVector<First_...> & vec, typename TupleVector<First_...>::DataType * pval_set)
      {
        static_assert(sizeof...(First_) == std::size_t(1), "invalid TupleVector size");
        _set_vec(vec.first(), pval_set);
      }
#else // all other compilers
      template<typename First_>
      static void _set_vec(const TupleVector<First_> & vec, typename TupleVector<First_>::DataType * pval_set)
      {
        _set_vec(vec.first(), pval_set);
      }
#endif

    public:
      /**
       * \brief Converts any vector to DenseVector
       *
       * \param[in] a The input vector.
       * \returns The vector a in DenseVector-format.
       */
      template<typename VT_>
      static DenseVector<typename VT_::MemType, typename VT_::DataType> value(const VT_ & a)
      {
        DenseVector<typename VT_::MemType, typename VT_::DataType> vec(a.size());
        auto * pvec(vec.elements());

        _set_vec(a, pvec);

        return vec;
      }

      /**
       * \brief Copy any vector into an allocated DenseVector
       *
       * \param[in] a The ouput vector of type DenseVector
       * \param[in] b The input vector
       */
      template<typename VT_>
      static void copy(DenseVector<typename VT_::MemType, typename VT_::DataType> & a, const VT_ & b)
      {
        if (a.size() != b.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vectors have not the same size!");

        auto * pa(a.elements());

        _set_vec(b, pa);
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_BLOCK_TO_SCALAR_HPP
