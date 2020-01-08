// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_PSD_HPP
#define KERNEL_SOLVER_PSD_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Steepest-Descent solver implementation
     *
     * This class implements a simple preconditioned steepest descent solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see Chapter 5.3.1, Algorithm 5.2 in \cite Saad03
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PSD :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<VectorType> BaseClass;

      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// temporary vectors
      VectorType _vec_r, _vec_q, _vec_z;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit PSD(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PSD", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] matrix
       * The system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] precond
       * The preconditioner. May be \c nullptr.
       */
      explicit PSD(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PSD", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      virtual String name() const override
      {
        return "PSD";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create three temporary vectors
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        _vec_z = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_z.clear();
        this->_vec_q.clear();
        this->_vec_r.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_r.copy(vec_def);

        // clear solution vector
        vec_cor.format();

        // apply solver
        this->_status = _apply_intern(vec_cor);

        // plot summary
        this->plot_summary();

        // return status
        return this->_status;
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply solver
        this->_status = _apply_intern(vec_sol);

        // plot summary
        this->plot_summary();

        // return status
        return this->_status;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol)
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_q(this->_vec_q);
        VectorType& vec_z(this->_vec_z);

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // apply preconditioner to defect vector
          // z[k] := M^{-1} * r[k]
          if(!this->_apply_precond(vec_z, vec_r, filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // q[k] := A*z[k]
          matrix.apply(vec_q, vec_z);
          filter.filter_def(vec_q);

          // alpha[k] := < z[k], r[k] > / < z[k], q[k] >
          DataType alpha = vec_z.dot(vec_r) / vec_z.dot(vec_q);

          // update solution vector:
          // x[k+1] := x[k] + alpha[k] * z[k]
          vec_sol.axpy(vec_z, vec_sol, alpha);

          // update defect vector:
          // r[k+1] := r[k] - alpha[k] * q[k]
          vec_r.axpy(vec_q, vec_r, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }
    }; // class PSD<...>

    /**
     * \brief Creates a new PSD solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new PSD object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new PSD solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new PSD object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PSD<Matrix_, Filter_>> new_psd(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PSD<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PSD_HPP
