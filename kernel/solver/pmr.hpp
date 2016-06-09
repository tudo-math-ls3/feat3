#pragma once
#ifndef KERNEL_SOLVER_PMR_HPP
#define KERNEL_SOLVER_PMR_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Minimal-Residual-Iteration solver implementation
     *
     * This class implements a simple preconditioned Minimal-Residual-Iteration solver.
     *
     * \attention
     * The <em>Minimal Residual <b>iteration</b></em> (MR) is \b not to be confused with the
     * <em>Minimal Residual <b>method</b></em> (MINRES), which is (mathematically) equivalent
     * to the Conjugate-Residual method!
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see Chapter 5.3.2, Algorithm 5.3 in \cite Saad03
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PMR :
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
      VectorType _vec_r, _vec_s, _vec_q, _vec_z;

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
      explicit PMR(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PMR", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      virtual String name() const override
      {
        return "PMR";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        _vec_r = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        _vec_z = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_z.clear();
        this->_vec_q.clear();
        this->_vec_s.clear();
        this->_vec_r.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_r.copy(vec_def);
        //this->_system_filter.filter_def(this->_vec_r);

        // clear solution vector
        vec_cor.format();

        // apply
        return _apply_intern(vec_cor, vec_def);
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_s(this->_vec_s);
        VectorType& vec_q(this->_vec_q);
        VectorType& vec_z(this->_vec_z);

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
          return status;

        // apply preconditioner to defect vector
        // s[0] := M^{-1} * r[0]
        if(!this->_apply_precond(vec_s, vec_r, filter))
          return Status::aborted;

        // start iterating
        while(status == Status::progress)
        {
          TimeStamp at;
          double mpi_start(Statistics::get_time_mpi_execute());

          // q[k] := A*s[k]
          matrix.apply(vec_q, vec_s);
          filter.filter_def(vec_q);

          // apply preconditioner to defect vector
          // z[k] := M^{-1} * q[k]
          if(!this->_apply_precond(vec_z, vec_q, filter))
            return Status::aborted;

          // alpha[k] := < q[k], s[k] > / < z[k], q[k] >
          DataType alpha = vec_q.dot(vec_s) / vec_z.dot(vec_q);

          // update solution vector:
          // x[k+1] := x[k] + alpha[k] * s[k]
          vec_sol.axpy(vec_s, vec_sol, alpha);

          // update defect vector:
          // r[k+1] := r[k] - alpha[k] * q[k]
          vec_r.axpy(vec_q, vec_r, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            double mpi_stop(Statistics::get_time_mpi_execute());
            Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
            return status;
          }

          // update preconditioned defect:
          // s[k+1] := s[k] - alpha[k] * z[k]
          vec_s.axpy(vec_z, vec_s, -alpha);

          TimeStamp bt;
          Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
          double mpi_stop(Statistics::get_time_mpi_execute());
          Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
        }

        // we should never reach this point...
        return Status::undefined;
      }
    }; // class PMR<...>

    /**
     * \brief Creates a new PMR solver object
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
     * A shared pointer to a new PMR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PMR<Matrix_, Filter_>> new_pmr(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PMR<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PMR<Matrix_, Filter_>> new_pmr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PMR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PMR<Matrix_, Filter_>> new_pmr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PMR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PMR_HPP
