#pragma once
#ifndef KERNEL_SOLVER_JACOBI_PRECOND_HPP
#define KERNEL_SOLVER_JACOBI_PRECOND_HPP 1

// includes, FEAT
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

namespace FEAT
{
  namespace Solver
  {
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
      /// Our base class
      typedef SolverBase<VectorType> BaseClass;

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
       * \param[in] filter
       * The system filter.
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

      /**
       * \brief Empty virtual destructor
       */
      virtual ~JacobiPrecond()
      {
      }

      /**
       * \brief Reads a solver configuration from a PropertyMap
       */
      virtual void read_config(PropertyMap* section) override
      {
        BaseClass::read_config(section);

        // Check if we have set _krylov_vim
        auto omega_p = section->query("omega");
        if(omega_p.second)
        {
          set_omega(DataType(std::stod(omega_p.first)));
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Jacobi";
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
        // extract matrix diagonal
        _matrix.extract_diag(_inv_diag);

        // invert diagonal elements
        _inv_diag.component_invert(_inv_diag, _omega);
      }

      /**
       * \brief Sets the damping parameter
       *
       * \param[in] omega
       * The new damping parameter.
       *
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        _omega = omega;
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
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_JACOBI_PRECOND_HPP
