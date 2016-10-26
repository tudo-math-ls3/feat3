#pragma once
#ifndef KERNEL_SOLVER_SSOR_PRECOND_HPP
#define KERNEL_SOLVER_SSOR_PRECOND_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>

namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      int cuda_ssor_forward_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      int cuda_ssor_backward_apply(int m, double * y, const double * x, double * csrVal, int * csrColInd, int ncolors, double omega,
          int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      void cuda_sor_done_symbolic(int * colored_row_ptr, int * rows_per_color, int * inverse_row_ptr);
      void cuda_sor_init_symbolic(int m, int nnz, double * csrVal, int * csrRowPtr, int * csrColInd, int & ncolors,
        int* & colored_row_ptr, int* & rows_per_color, int* & inverse_row_ptr);
    }

    /**
     * \brief SSOR preconditioner implementation
     *
     * This class implements a simple SSOR preconditioner,
     * e.g. zero fill-in and no pivoting.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - LAFEM::SparseMatrixCSR
     * - LAFEM::SparseMatrixELL
     *
     * Moreover, this implementation supports only Mem::Main
     *
     * \author Dirk Ribbrock
     */
    template<typename Matrix_, typename Filter_>
    class SSORPrecond :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef typename MatrixType::IndexType IndexType;

    protected:
      const MatrixType& _matrix;
      const FilterType& _filter;
      DataType _omega;

      void _apply_intern(const LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType>& matrix, VectorType& vec_cor, const VectorType& vec_def)
      {
        // create pointers
        DataType * pout(vec_cor.elements());
        const DataType * pin(vec_def.elements());
        const DataType * pval(matrix.val());
        const IndexType * pcol_ind(matrix.col_ind());
        const IndexType * prow_ptr(matrix.row_ptr());
        const IndexType n((IndexType(matrix.rows())));

        // __forward-insertion__
        // iteration over all rows
        for (Index i(0); i < n; ++i)
        {
          IndexType col;
          DataType d(0);
          // iteration over all elements on the left side of the main-diagonal
          for (col = prow_ptr[i]; pcol_ind[col] < i; ++col)
          {
            d += pval[col] * pout[pcol_ind[col]];
          }
          pout[i] = (pin[i] - _omega * d) / pval[col];
        }

        // __backward-insertion__
        // iteration over all rows
        for (Index i(n); i > 0;)
        {
          --i;
          IndexType col;
          DataType d(0);
          // iteration over all elements on the right side of the main-diagonal
          for (col = prow_ptr[i+1] - IndexType(1); pcol_ind[col] > i; --col)
          {
            d += pval[col] * pout[pcol_ind[col]];
          }
          pout[i] -= _omega * d / pval[col];
        }
      }

      void _apply_intern(const LAFEM::SparseMatrixELL<Mem::Main, DataType, IndexType>& matrix, VectorType& vec_cor, const VectorType& vec_def)
      {
        // create pointers
        DataType * pout(vec_cor.elements());
        const DataType * pin(vec_def.elements());
        const DataType * pval(matrix.val());
        const IndexType * pcol_ind(matrix.col_ind());
        const IndexType * pcs(matrix.cs());
        const IndexType * prl(matrix.rl());
        const IndexType C((IndexType(matrix.C())));
        const IndexType n((IndexType(matrix.rows())));

        // __forward-insertion__
        // iteration over all rows
        for (IndexType i(0); i < n; ++i)
        {
          IndexType col;
          DataType d(0);
          // iteration over all elements on the left side of the main-diagonal
          for (col = pcs[i/C] + i%C; pcol_ind[col] < i; col += C)
          {
            d += pval[col] * pout[pcol_ind[col]];
          }
          pout[i] = (pin[i] - _omega * d) / pval[col];
        }

        // __backward-insertion__
        // iteration over all rows
        for (IndexType i(n); i > 0;)
        {
          --i;
          IndexType col;
          DataType d(0);
          // iteration over all elements on the right side of the main-diagonal
          for (col = pcs[i/C] + i%C + C * (prl[i] - 1); pcol_ind[col] > i; col -= C)
          {
            d += pval[col] * pout[pcol_ind[col]];
          }
          pout[i] -= _omega * d / pval[col];
        }
      }

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The filter to be used for the correction vector.
       *
       * \param[in] omega
       * Damping
       *
       */
      explicit SSORPrecond(const MatrixType& matrix, const FilterType& filter, const DataType omega = DataType(1)) :
        _matrix(matrix),
        _filter(filter),
        _omega(omega)
      {
        if (_matrix.columns() != _matrix.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "SSOR";
      }

      virtual void init_symbolic() override
      {
      }

      virtual void done_symbolic() override
      {
      }

      virtual void init_numeric() override
      {
      }

      virtual void done_numeric() override
      {
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        XASSERTM(_matrix.rows() == vec_cor.size(), "matrix / vector size mismatch!");
        XASSERTM(_matrix.rows() == vec_def.size(), "matrix / vector size mismatch!");

        TimeStamp ts_start;

        // copy in-vector to out-vector
        vec_cor.copy(vec_def);

        _apply_intern(_matrix, vec_cor, vec_def);

        vec_cor.scale(vec_cor, _omega * (DataType(2.0) - _omega));

        this->_filter.filter_cor(vec_cor);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_matrix.used_elements() *2 + 6 * vec_cor.size());

        return Status::success;
      }
    }; // class SSORPrecond<SparseMatrixCSR<Mem::Main>>

    /**
     * \brief SSOR preconditioner implementation
     *
     * This class implements a simple SSOR preconditioner,
     * e.g. zero fill-in and no pivoting.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - LAFEM::SparseMatrixCSR
     *
     * Moreover, this implementation supports only CUDA and uint containers.
     *
     * \note This class need at least cuda version 7.
     *
     * \author Dirk Ribbrock
     */
    template<typename Filter_>
    class SSORPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int>::VectorTypeL>
    {
    public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned int> MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;

    protected:
      const MatrixType& _matrix;
      const FilterType& _filter;
      double _omega;
      // row ptr permutation, sorted by color(each color sorted by amount of rows), start/end index per row
      int * _colored_row_ptr;
      // amount of rows per color (sorted by amount of rows)
      int * _rows_per_color;
      // mapping of idx to native row number
      int * _inverse_row_ptr;
      // number of colors
      int _ncolors;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The filter to be used for the correction vector.
       *
       * \param[in] omega
       * Damping
       *
       */
      explicit SSORPrecond(const MatrixType& matrix, const FilterType& filter, const DataType omega = DataType(1)) :
        _matrix(matrix),
        _filter(filter),
        _omega(omega),
        _colored_row_ptr(nullptr),
        _rows_per_color(nullptr),
        _inverse_row_ptr(nullptr),
        _ncolors(0)
      {
        if (_matrix.columns() != _matrix.rows())
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Matrix is not square!");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "SSOR";
      }

      virtual void init_symbolic() override
      {
        Intern::cuda_sor_init_symbolic((int)_matrix.rows(), (int)_matrix.used_elements(), (double*)_matrix.val(), (int*)_matrix.row_ptr(), (int*)_matrix.col_ind(), _ncolors,
          _colored_row_ptr, _rows_per_color, _inverse_row_ptr);
      }

      virtual void done_symbolic() override
      {
        Intern::cuda_sor_done_symbolic(_colored_row_ptr, _rows_per_color, _inverse_row_ptr);
      }

      virtual void init_numeric() override
      {
      }

      virtual void done_numeric() override
      {
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        XASSERTM(_matrix.rows() == vec_cor.size(), "matrix / vector size mismatch!");
        XASSERTM(_matrix.rows() == vec_def.size(), "matrix / vector size mismatch!");

        TimeStamp ts_start;

        int status = Intern::cuda_ssor_forward_apply((int)vec_cor.size(), vec_cor.elements(), vec_def.elements(), (double*)_matrix.val(), (int*)_matrix.col_ind(), _ncolors, _omega, _colored_row_ptr, _rows_per_color, _inverse_row_ptr);
        status |= Intern::cuda_ssor_backward_apply((int)vec_cor.size(), vec_cor.elements(), vec_def.elements(), (double*)_matrix.val(), (int*)_matrix.col_ind(), _ncolors, _omega, _colored_row_ptr, _rows_per_color, _inverse_row_ptr);

        vec_cor.scale(vec_cor, _omega * (2.0 - _omega));

        this->_filter.filter_cor(vec_cor);

        TimeStamp ts_stop;
        Statistics::add_time_precon(ts_stop.elapsed(ts_start));
        Statistics::add_flops(_matrix.used_elements() *2 + 6 * vec_cor.size());

        return (status == 0) ? Status::success :  Status::aborted;
      }
    }; // class SSORPrecond<SparseMatrixCSR<Mem::CUDA>>

    template<typename Filter_>
    class SSORPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned int> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::DataType DataType;

      explicit SSORPrecond(const MatrixType&, const FilterType&, const DataType = DataType(1))
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }

      String name() const override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }
    };

    template<typename Filter_>
    class SSORPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, float, unsigned long> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::DataType DataType;

      explicit SSORPrecond(const MatrixType&, const FilterType&, const DataType = DataType(1))
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }

      String name() const override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }
    };

    template<typename Filter_>
    class SSORPrecond<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long>, Filter_> :
      public SolverBase<LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long>::VectorTypeL>
    {
      public:
      typedef LAFEM::SparseMatrixCSR<Mem::CUDA, double, unsigned long> MatrixType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::DataType DataType;

      explicit SSORPrecond(const MatrixType&, const FilterType&, const DataType = DataType(1))
      {
      }

      Status apply(VectorType &, const VectorType &) override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }

      String name() const override
      {
          throw InternalError(__func__, __FILE__, __LINE__, "not implemented yet!");
      }
    };

    /**
     * \brief Creates a new SSORPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The damping value.
     *
     * \returns
     * A shared pointer to a new SSORPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<SSORPrecond<Matrix_, Filter_>> new_ssor_precond(
      const Matrix_& matrix, const Filter_& filter, const typename Matrix_::DataType omega = 1.0)
    {
      return std::make_shared<SSORPrecond<Matrix_, Filter_>>
        (matrix, filter, omega);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SSOR_PRECOND_HPP
