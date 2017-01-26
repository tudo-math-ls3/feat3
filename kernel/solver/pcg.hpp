#pragma once
#ifndef KERNEL_SOLVER_PCG_HPP
#define KERNEL_SOLVER_PCG_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Conjugate-Gradient solver implementation
     *
     * This class implements a simple PCG solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see Chapter 9.2, Algorithm 9.1 in \cite Saad03
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class PCG :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      /// The type of matrix this solver can be applied to
      typedef Matrix_ MatrixType;
      /// The filter for projecting solution, rhs, defect and correction vectors to subspaces
      typedef Filter_ FilterType;
      /// The vector type this solver can be applied to
      typedef typename MatrixType::VectorTypeR VectorType;
      /// The floating point precision
      typedef typename MatrixType::DataType DataType;
      /// Our base class
      typedef PreconditionedIterativeSolver<VectorType> BaseClass;
      /// The type of the preconditioner that can be used
      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// The defect vector
      VectorType _vec_r;
      /// The update (or search) direction
      VectorType _vec_p;
      /// Temporary vector, used for e.g. the preconditioned defect
      VectorType _vec_t;

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
      explicit PCG(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PCG", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
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
     *
     */
      explicit PCG(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PCG", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      /// \copydoc SolverBase::name()
      virtual String name() const override
      {
        return "PCG";
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create three temporary vectors
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_p = this->_system_matrix.create_vector_r();
        _vec_t = this->_system_matrix.create_vector_r();
      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        this->_vec_t.clear();
        this->_vec_p.clear();
        this->_vec_r.clear();
        BaseClass::done_symbolic();
      }

      /// \copydoc SolverBase::apply()
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

      /// \copydoc IterativeSolver::correct()
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      /**
       * \brief Internal function, applies the solver
       *
       * \param[in] vec_sol
       * The current solution vector, gets overwritten
       *
       * \param[in] vec_rhs
       * The right hand side vector. This is unused in this function, as the intial defect was alredy computed and
       * stored in _vec_r.
       *
       * \returns A status code.
       */
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_p(this->_vec_p);
        // Note: q and z are temporary vectors whose usage does
        // not overlap, so we use the same vector for them
        VectorType& vec_q(this->_vec_t);
        VectorType& vec_z(this->_vec_t);

        // Note:
        // In the algorithm below, the following relations hold:
        // q[k] = A * p[k]
        // z[k] = M^{-1} * r[k]

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        // apply preconditioner to defect vector
        // p[0] := M^{-1} * r[0]
        if(!this->_apply_precond(vec_p, vec_r, filter))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        // compute initial gamma:
        // gamma[0] := < r[0], p[0] >
        DataType gamma = vec_r.dot(vec_p);

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // q[k] := A*p[k]
          matrix.apply(vec_q, vec_p);
          filter.filter_def(vec_q);

          // compute alpha
          // alpha[k] := gamma[k] / < q[k], p[k] >
          DataType alpha = gamma / vec_q.dot(vec_p);

          // update solution vector:
          // x[k+1] := x[k] + alpha[k] * p[k]
          vec_sol.axpy(vec_p, vec_sol, alpha);

          // update defect vector:
          // r[k+1] := r[k] - alpha[k] * q[k]
          vec_r.axpy(vec_q, vec_r, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // apply preconditioner
          // z[k+1] := M^{-1} * r[k+1]
          if(!this->_apply_precond(vec_z, vec_r, filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // compute new gamma:
          // gamma[k+1] := < r[k+1] , z[k+1] >
          DataType gamma2 = gamma;
          gamma = vec_r.dot(vec_z);

          // compute beta:
          // beta[k] := gamma[k+1] / gamma[k]
          DataType beta = gamma / gamma2;

          // update direction vector:
          // p[k+1] := z[k+1] + beta[k] * p[k]
          vec_p.axpy(vec_p, vec_z, beta);
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }
    }; // class PCG<...>

    /**
     * \brief Creates a new PCG solver object
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
     * A shared pointer to a new PCG object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new PCG solver object using a PropertyMap
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
     * A shared pointer to a new PCG object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCG<Matrix_, Filter_>> new_pcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PCG<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PCG_HPP
