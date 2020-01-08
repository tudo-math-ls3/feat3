// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_PIPEPCG_HPP
#define KERNEL_SOLVER_PIPEPCG_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) pipelined Conjugate-Gradient solver implementation from Ghysels and Vnaroose
     *
     * This method has only a single non-blocking reduction per iteration, compared to 2 blocking for standard CG.  The
     * non-blocking reduction is overlapped by the matrix-vector product and preconditioner application.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * Reference:
     * P. Ghysels and W. Vanroose, "Hiding global synchronization latency in the preconditioned Conjugate Gradient algorithm",
     *
     * \see https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/impls/cg/pipecg/pipecg.c.html#KSPPIPECG
     *
     *
     * \author Dirk Ribbrock
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PipePCG :
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
      VectorType _vec_r, _vec_u, _vec_w, _vec_z, _vec_q, _vec_s, _vec_p, _vec_m, _vec_n;

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
      explicit PipePCG(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PipePCG", precond),
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
      explicit PipePCG(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("PipePCG", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      virtual String name() const override
      {
        return "PipePCG";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create temporary vectors
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_u = this->_system_matrix.create_vector_r();
        _vec_w = this->_system_matrix.create_vector_r();
        _vec_z = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_p = this->_system_matrix.create_vector_r();
        _vec_m = this->_system_matrix.create_vector_r();
        _vec_n = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_r.clear();
        this->_vec_u.clear();
        this->_vec_w.clear();
        this->_vec_z.clear();
        this->_vec_q.clear();
        this->_vec_s.clear();
        this->_vec_p.clear();
        this->_vec_m.clear();
        this->_vec_n.clear();
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
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_u(this->_vec_u);
        VectorType& vec_w(this->_vec_w);
        VectorType& vec_z(this->_vec_z);
        VectorType& vec_q(this->_vec_q);
        VectorType& vec_s(this->_vec_s);
        VectorType& vec_p(this->_vec_p);
        VectorType& vec_m(this->_vec_m);
        VectorType& vec_n(this->_vec_n);

        DataType gamma, gamma_old, delta, beta, alpha;
        gamma = DataType(0);
        gamma_old = DataType(0);
        delta = DataType(0);
        beta = DataType(0);
        alpha = DataType(0);

        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        if(!this->_apply_precond(vec_u, vec_r, filter))
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        matrix.apply(vec_w, vec_u);
        filter.filter_def(vec_w);

        pre_iter.destroy();

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          auto norm_def_cur = vec_r.norm2_async();
          auto dot_gamma = vec_r.dot_async(vec_u);
          auto dot_delta = vec_w.dot_async(vec_u);

          if(!this->_apply_precond(vec_m, vec_w, filter))
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          matrix.apply(vec_n, vec_m);
          filter.filter_def(vec_n);

          gamma = dot_gamma->wait();
          delta = dot_delta->wait();

          status = this->_update_defect(norm_def_cur->wait());
          if(status != Status::progress)
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          if (this->_num_iter == 1) // num_iter has already been increased to 1 by previous _update_defect call
          {
            alpha = gamma / delta;
            vec_z.copy(vec_n);
            vec_q.copy(vec_m);
            vec_p.copy(vec_u);
            vec_s.copy(vec_w);
          }
          else
          {
            beta = gamma / gamma_old;
            alpha = gamma / (delta - beta / alpha * gamma);
            vec_z.axpy(vec_z, vec_n, beta);
            vec_q.axpy(vec_q, vec_m, beta);
            vec_p.axpy(vec_p, vec_u, beta);
            vec_s.axpy(vec_s, vec_w, beta);
          }

          vec_sol.axpy(vec_p, vec_sol, alpha);

          /* Stabilization, recalculate real values
          if(this->_num_iter > 40 && this->_num_iter % 50 == 0)
          {
            matrix.apply(vec_r, vec_sol, DataType(-1));
            filter.filter_def(vec_r);
            if(!this->_apply_precond(vec_u, vec_r, filter))
              return Status::aborted;
            matrix.apply(vec_w, vec_u);
            filter.filter_def(vec_w);
          }
          else */
          {
            vec_u.axpy(vec_q, vec_u, -alpha);
            vec_w.axpy(vec_z, vec_w, -alpha);
            vec_r.axpy(vec_s, vec_r, -alpha);
          }

          gamma_old = gamma;
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }

    }; // class PipePCG<...>

    /**
     * \brief Creates a new PipePCG solver object
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
     * A shared pointer to a new PipePCG object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new PipePCG solver object using a PropertyMap
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
     * A shared pointer to a new PipePCG object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter, std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PipePCG<Matrix_, Filter_>> new_pipepcg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<PipePCG<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PIPEPCG_HPP
