#pragma once
#ifndef KERNEL_SOLVER_JACOBI_PRECOND_HPP
#define KERNEL_SOLVER_JACOBI_PRECOND_HPP 1

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
    /// \cond internal
    namespace Intern
    {
      template<typename MT_>
      struct JacobiHelper;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Jacobi preconditioner implementation
     *
     * This class implements a simple damped Jacobi preconditioner.
     *
     * This implementation works for the following matrix types and combinations thereof:
     * - all LAFEM::SparseMatrix types
     * - LAFEM::DenseMatrix
     * - LAFEM::PowerDiagMatrix
     * - LAFEM::PowerFullMatrix
     * - LAFEM::TupleDiagMatrix
     * - Global::Matrix
     *
     * Moreover, this implementation supports all Mem architectures, as well as all
     * data and index types.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class JacobiPrecond :
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
      DataType _omega;
      VectorType _inv_diag;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The matrix whose main diagonal is to be used.
       *
       * \param[in] omega
       * The damping parameter for the preconditioner.
       */
      explicit JacobiPrecond(const MatrixType& matrix, const FilterType& filter, DataType omega = DataType(1)) :
        _matrix(matrix),
        _filter(filter),
        _omega(omega)
      {
      }

      virtual void init_symbolic() override
      {
        _inv_diag = _matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        _inv_diag.clear();
      }

      virtual void init_numeric() override
      {
        // Note: extraction and inversion of the main diagonal is split up into
        // two steps, as it would not be possible to support Global::Matrix otherwise.

        // extract matrix diagonal
        Intern::JacobiHelper<MatrixType>::extract_diag(_inv_diag, _matrix);

        // invert diagonal elements
        _inv_diag.component_invert(_inv_diag, _omega);
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        vec_cor.component_product(_inv_diag, vec_def);
        this->_filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class JacobiPrecond<...>

    /**
     * \brief Creates a new JacobiPrecond solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The damping parameter for Jacobi.
     *
     * \returns
     * A shared pointer to a new JacobiPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<JacobiPrecond<Matrix_, Filter_>> new_jacobi_precond(
      const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1))
    {
      return std::make_shared<JacobiPrecond<Matrix_, Filter_>>(matrix, filter, omega);
    }

    /// \cond internal
    namespace Intern
    {
      // generic JacobiHelper implemenation for scalar matrices
      template<typename MT_>
      struct JacobiHelper
      {
        static void extract_diag(typename MT_::VectorTypeL& diag, const MT_& matrix)
        {
          matrix.extract_diag(diag);
        }
      };

      // JacobiHelper specialisation for PowerDiagMatrix
      template<typename SMT_, int n_>
      struct JacobiHelper< LAFEM::PowerDiagMatrix<SMT_, n_> >
      {
        static_assert(n_ > 1, "invalid block size");
        typedef LAFEM::PowerDiagMatrix<SMT_, n_> MT;
        static void extract_diag(typename MT::VectorTypeL& diag, const MT& matrix)
        {
          JacobiHelper<SMT_>::extract_diag(diag.first(), matrix.first());
          JacobiHelper<LAFEM::PowerDiagMatrix<SMT_, n_-1>>::extract_diag(diag.rest(), matrix.rest());
        }
      };

      template<typename SMT_>
      struct JacobiHelper< LAFEM::PowerDiagMatrix<SMT_, 1> >
      {
        typedef LAFEM::PowerDiagMatrix<SMT_, 1> MT;
        static void extract_diag(typename MT::VectorTypeL& diag, const MT& matrix)
        {
          JacobiHelper<SMT_>::extract_diag(diag.first(), matrix.first());
        }
      };

      // JacobiHelper specialisation for PowerFullMatrix
      template<typename SMT_, int n_>
      struct JacobiHelper< LAFEM::PowerFullMatrix<SMT_, n_, n_> >
      {
        static_assert(n_ > 1, "invalid block size");
        template<typename MT_>
        static void extract_diag(typename MT_::VectorTypeL& diag, const MT_& matrix)
        {
          JacobiHelper<LAFEM::PowerFullMatrix<SMT_, n_-1, n_-1>>::extract_diag(diag, matrix);
          JacobiHelper<SMT_>::extract_diag(diag.template at<n_-1>(), matrix.template at<n_-1,n_-1>());
        }
      };

      template<typename SMT_>
      struct JacobiHelper< LAFEM::PowerFullMatrix<SMT_, 1, 1> >
      {
        template<typename MT_>
        static void extract_diag(typename MT_::VectorTypeL& diag, const MT_& matrix)
        {
          JacobiHelper<SMT_>::extract_diag(diag.template at<0>(), matrix.template at<0,0>());
        }
      };

      // JacobiHelper specialisation for TupleDiagMatrix
      template<typename First_, typename... Rest_>
      struct JacobiHelper< LAFEM::TupleDiagMatrix<First_, Rest_...> >
      {
        typedef LAFEM::TupleDiagMatrix<First_, Rest_...> MT;
        static void extract_diag(typename MT::VectorTypeL& diag, const MT& matrix)
        {
          JacobiHelper<First_>::extract_diag(diag.first(), matrix.first());
          JacobiHelper<LAFEM::TupleDiagMatrix<Rest_...>>::extract_diag(diag.rest(), matrix.rest());
        }
      };

      template<typename First_>
      struct JacobiHelper< LAFEM::TupleDiagMatrix<First_> >
      {
        typedef LAFEM::TupleDiagMatrix<First_> MT;
        static void extract_diag(typename MT::VectorTypeL& diag, const MT& matrix)
        {
          JacobiHelper<First_>::extract_diag(diag.first(), matrix.first());
        }
      };

      // JacobiHelper specialisation for Global::Matrix
      template<typename LocalMatrix_>
      struct JacobiHelper< Global::Matrix<LocalMatrix_> >
      {
        typedef Global::Matrix<LocalMatrix_> MT;
        static void extract_diag(typename MT::VectorTypeL& diag, const MT& matrix)
        {
          // extract local diagonal
          JacobiHelper<LocalMatrix_>::extract_diag(*diag, *matrix);
          // synchronise to convert from type-0 to type-1
          diag.sync_0();
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_JACOBI_PRECOND_HPP
