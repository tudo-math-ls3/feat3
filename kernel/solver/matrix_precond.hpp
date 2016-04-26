#pragma once
#ifndef KERNEL_SOLVER_MATRIX_PRECOND_HPP
#define KERNEL_SOLVER_MATRIX_PRECOND_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/tuple_diag_matrix.hpp>
#include <kernel/global/matrix.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Matrix preconditioner implementation
     *
     * This class wrapps a matrix to be used as a preconditioner
     *
     * \author Dirk Ribbrock
     */
    template<typename Matrix_, typename Filter_>
    class MatrixPrecond :
      public SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeL VectorType;
      typedef typename MatrixType::DataType DataType;

    protected:
      const MatrixType& _matrix;
      const FilterType& _filter;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix to be used.
       *
       * \param[in] filter
       * The system filter.
       */
      explicit MatrixPrecond(const MatrixType& matrix, const FilterType& filter) :
        _matrix(matrix),
        _filter(filter)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "MatrixPrecond";
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

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        _matrix.apply(vec_cor, vec_def);
        this->_filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class MatrixPrecond<...>

    /**
     * \brief Creates a new MatrixPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new JacobiPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<MatrixPrecond<Matrix_, Filter_>> new_matrix_precond(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<MatrixPrecond<Matrix_, Filter_>>(matrix, filter);
    }
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_MATRIX_PRECOND_HPP
