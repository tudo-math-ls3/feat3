#pragma once
#ifndef KERNEL_SOLVER_PCR_HPP
#define KERNEL_SOLVER_PCR_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Conjugate-Residual solver implementation
     *
     * This implementation is based on the unpreconditioned algorithm from \cite Saad03.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see Chapter 6.8, Algorithm 6.20 in \cite Saad03
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PCR :
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
      /// vectors
      VectorType _vec_r, _vec_s, _vec_p, _vec_q, _vec_t;

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
      explicit PCR(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PCR", precond),
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
       *
       */
      explicit PCR(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PCR", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      virtual String name() const override
      {
        return "PCR";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create five temporary vectors
        _vec_p = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_t = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_t.clear();
        this->_vec_s.clear();
        this->_vec_r.clear();
        this->_vec_q.clear();
        this->_vec_p.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_r.copy(vec_def);

        // clear solution vector
        vec_cor.format();

        // apply
        Status st(_apply_intern(vec_cor, vec_def));
        this->plot_summary(st);
        return st;
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply
        Status st(_apply_intern(vec_sol, vec_rhs));
        this->plot_summary(st);
        return st;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_p(this->_vec_p);
        VectorType& vec_q(this->_vec_q);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_s(this->_vec_s);
        // Note: y and z are temporary vectors whose usage does
        // not overlap, so we use the same vector for them
        VectorType& vec_y(this->_vec_t);
        VectorType& vec_z(this->_vec_t);

        // Note:
        // In the algorithm below, the following relations hold:
        // q[k] = A * p[k]
        // s[k] = M^{-1} * r[k]
        // y[k] = A * s[k]
        // z[k] = M^{-1} * q[k]

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        // apply preconditioner to defect vector:
        // s[0] := M^{-1} * r[0]
        if(!this->_apply_precond(vec_s, vec_r, filter))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        // compute initial p and q
        // p[0] := s[0]
        vec_p.copy(vec_s);

        // q[0] := A*p[0]
        matrix.apply(vec_q, vec_p);
        filter.filter_def(vec_q);

        // compute initial gamma:
        // gamma[0] := < p[0], q[0] >
        DataType gamma = vec_p.dot(vec_q);

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // z[k] := M^{-1} * q[k]
          if(!this->_apply_precond(vec_z, vec_q, filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // compute alpha
          // alpha[k] := gamma[k] / < z[k], q[k] >
          DataType alpha = gamma / vec_z.dot(vec_q);

          // update solution vector
          // x[k+1] := x[k] + alpha[k] * p[k]
          vec_sol.axpy(vec_p, vec_sol, alpha);

          // update real residual vector
          // r[k+1] := r[k] - alpha[k] * q[k]
          vec_r.axpy(vec_q, vec_r, -alpha);

          // compute defect norm and check for convergence
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // update preconditioned residual vector
          // s[k+1] := s[k] - alpha[k] * z[k]
          vec_s.axpy(vec_z, vec_s, -alpha);

          // y[k+1] := A*s[k+1]
          matrix.apply(vec_y, vec_s);
          filter.filter_def(vec_y);

          // compute new gamma:
          // gamma[k+1] := < s[k+1], y[k+1] >
          DataType gamma2 = gamma;
          gamma = vec_s.dot(vec_y);

          // compute beta:
          // beta[k] := gamma[k+1] / gamma[k]
          DataType beta = gamma / gamma2;

          // update direction vectors:
          // p[k+1] := s[k+1] + beta[k] * p[k]
          vec_p.axpy(vec_p, vec_s, beta);

          // q[k+1] := y[k+1] + beta[k] * q[k]
          vec_q.axpy(vec_q, vec_y, beta);
        }

        // we should never reach this point...
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }
      }
    }; // class PCR<...>

    /**
     * \brief Creates a new PCR solver object
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
     * A shared pointer to a new PCR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new PCR solver object using a PropertyMap
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
     * A shared pointer to a new PCR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCR<Matrix_, Filter_>> new_pcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PCR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PCR_HPP
