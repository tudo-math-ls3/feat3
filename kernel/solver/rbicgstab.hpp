// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_RBICGSTAB_HPP
#define KERNEL_SOLVER_RBICGSTAB_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) reordered BiCGStab solver implementation
     *
     * This class implements a reordered BiCGStab solver, which hides
     * the global synchronisation behind preconditioner calls.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see \cite Krasnopolsky10
     *
     * \author Dirk Ribbrock
     */
    template<typename Matrix_, typename Filter_>
    class RBiCGStab :
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
      /// temp vectors
      VectorType _vec_v, _vec_vh, _vec_z, _vec_s, _vec_t, _vec_th, _vec_r, _vec_r0;

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
      explicit RBiCGStab(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("RBiCGStab", precond),
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
      explicit RBiCGStab(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("RBiCGStab", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      /// \copydoc SolverBase::name()
      virtual String name() const override
      {
        return "RBiCGStab";
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create fife temporary vectors
        _vec_v = this->_system_matrix.create_vector_r();
        _vec_vh = this->_system_matrix.create_vector_r();
        _vec_z = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_t = this->_system_matrix.create_vector_r();
        _vec_th = this->_system_matrix.create_vector_r();
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_r0 = this->_system_matrix.create_vector_r();

      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        this->_vec_v.clear();
        this->_vec_vh.clear();
        this->_vec_z.clear();
        this->_vec_s.clear();
        this->_vec_t.clear();
        this->_vec_th.clear();
        this->_vec_r.clear();
        this->_vec_r0.clear();
        BaseClass::done_symbolic();
      }

      /// \copydoc SolverBase::apply()
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

      /// \copydoc IterativeSolver::correct()
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
      /**
       * \brief Internal function, applies the solver
       *
       * \param[in] vec_sol
       * The current solution vector, gets overwritten
       *
       * \param[in] vec_rhs
       * The right hand side vector. This is unused in this function, as the initial defect was alredy computed and
       * stored in _vec_r.
       *
       * \returns A status code.
       */
      virtual Status _apply_intern(VectorType& vec_sol)
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_v(this->_vec_v);
        VectorType& vec_vh(this->_vec_vh);
        VectorType& vec_z(this->_vec_z);
        VectorType& vec_s(this->_vec_s);
        VectorType& vec_t(this->_vec_t);
        VectorType& vec_th(this->_vec_th);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_r0(this->_vec_r0);

        DataType rho(0);
        DataType rho_old(0);
        DataType delta(0);
        DataType alpha(0);
        DataType theta(0);
        DataType omega(0);
        DataType beta(0);
        DataType phi(0);
        DataType psi(0);

        vec_r0.copy(vec_r);

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        auto dot_rho = vec_r.dot_async(vec_r);

        // apply preconditioner to defect vector
        if(!this->_apply_precond(vec_z, vec_r, filter))
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        vec_vh.copy(vec_z);

        rho = dot_rho->wait();
        rho_old = rho;

        pre_iter.destroy();

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          matrix.apply(vec_v, vec_vh);
          filter.filter_def(vec_v);

          auto dot_delta = vec_v.dot_async(vec_r0);

          // apply preconditioner
          if(!this->_apply_precond(vec_s, vec_v, filter))
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          delta = dot_delta->wait();
          alpha = rho / delta;

          vec_th.axpy(vec_s, vec_z, -alpha);

          matrix.apply(vec_t, vec_th);
          filter.filter_def(vec_t);

          vec_sol.axpy(vec_vh, vec_sol, alpha);
          vec_r.axpy(vec_v, vec_r, -alpha);
          auto norm_def_half = vec_r.norm2_async();

          auto dot_theta = vec_t.dot_async(vec_r);
          auto dot_phi = vec_t.dot_async(vec_t);
          auto dot_psi = vec_t.dot_async(vec_r0);

          // apply preconditioner
          if(!this->_apply_precond(vec_z, vec_t, filter))
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          theta = dot_theta->wait();
          phi = dot_phi->wait();
          psi = dot_psi->wait();

          // Check if we are already converged or failed after the "half" update
          {
            Status status_half(Status::progress);

            DataType def_old(this->_def_cur);
            DataType def_half(norm_def_half->wait());

            // ensure that the defect is neither NaN nor infinity
            if(!Math::isfinite(def_half))
            {
              status_half = Status::aborted;
            }
            // is diverged?
            else if(this->is_diverged(def_half))
            {
              this->_def_cur = def_half;
              status_half = Status::diverged;
            }
            // is converged?
            else if(this->is_converged(def_half))
            {
              this->_def_cur = def_half;
              status_half = Status::success;
            }

            // If we break early, we still count this iteration and plot it if necessary
            if(status_half != Status::progress)
            {
              ++this->_num_iter;
              this->_def_prev = def_old;
              // plot?
              if(this->_plot_iter())
                this->_plot_iter_line(this->_num_iter, this->_def_cur, def_old);

              stat.destroy();
              Statistics::add_solver_expression(
                  std::make_shared<ExpressionEndSolve>(this->name(), status_half, this->get_num_iter()));

              return status_half;
            }
          }

          omega = theta / phi;

          vec_sol.axpy(vec_th, vec_sol, omega);
          vec_r.axpy(vec_t, vec_r, -omega);

          auto norm_def_cur = vec_r.norm2_async();

          rho = -omega * psi;
          vec_z.axpy(vec_z, vec_th, -omega);
          beta = (rho / rho_old) * (alpha / omega);
          rho_old = rho;
          vec_vh.axpy(vec_s, vec_vh, -omega);
          vec_vh.axpy(vec_vh, vec_z, beta);

          status = this->_update_defect(norm_def_cur->wait());
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
    }; // class RBiCGStab<...>

    /**
     * \brief Creates a new RBiCGStab solver object
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
     * A shared pointer to a new RBiCGStab object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new RBiCGStab solver object using a PropertyMap
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
     * A shared pointer to a new RBiCGStab object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RBiCGStab<Matrix_, Filter_>> new_rbicgstab(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<RBiCGStab<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_RBICGSTAB_HPP
