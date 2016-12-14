#pragma once
#ifndef KERNEL_SOLVER_GROPPPCG_HPP
#define KERNEL_SOLVER_GROPPPCG_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) pipelined Conjugate-Gradient solver implementation from Bill Gropp
     *
     * This pcg method has two reductions, one of which is overlapped with the matrix-vector product and one of which is
     * overlapped with the preconditioner.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see http://www.cs.uiuc.edu/~wgropp/bib/talks/tdata/2012/icerm.pdf
     * \see https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/impls/cg/groppcg/groppcg.c.html#KSPGROPPCG
     *
     * \author Dirk Ribbrock
     */
    template<
      typename Matrix_,
      typename Filter_>
    class GroppPCG :
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
      VectorType _vec_r, _vec_z, _vec_s, _vec_p, _vec_S, _vec_Z;

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
      explicit GroppPCG(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("GroppPCG", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      explicit GroppPCG(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("GroppPCG", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      virtual String name() const override
      {
        return "GroppPCG";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create three temporary vectors
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_z = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_p = this->_system_matrix.create_vector_r();
        _vec_S = this->_system_matrix.create_vector_r();
        _vec_Z = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_r.clear();
        this->_vec_z.clear();
        this->_vec_s.clear();
        this->_vec_p.clear();
        this->_vec_S.clear();
        this->_vec_Z.clear();
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
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_z(this->_vec_z);
        VectorType& vec_s(this->_vec_s);
        VectorType& vec_p(this->_vec_p);
        VectorType& vec_S(this->_vec_S);
        VectorType& vec_Z(this->_vec_Z);

        DataType gamma, gamma_new, beta, alpha, t;
        gamma = DataType(0);
        gamma_new = DataType(0);
        beta = DataType(0);
        alpha = DataType(0);
        t = DataType(0);

        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        if(!this->_apply_precond(vec_z, vec_r, filter))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        vec_p.copy(vec_z);

        gamma = vec_z.dot(vec_r);

        matrix.apply(vec_s, vec_p);
        filter.filter_def(vec_s);

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          auto dot_t = vec_s.dot_async(vec_p);
          if(!this->_apply_precond(vec_S, vec_s, filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }
          t = dot_t->wait();

          alpha = gamma / t;

          vec_r.axpy(vec_s, vec_r, -alpha);
          auto norm_def_cur = vec_r.norm2_async();

          vec_sol.axpy(vec_p, vec_sol, alpha);

          vec_z.axpy(vec_S, vec_z, -alpha);

          auto dot_gamma_new = vec_r.dot_async(vec_z);
          matrix.apply(vec_Z, vec_z);
          filter.filter_def(vec_Z);

          gamma_new = dot_gamma_new->wait();
          beta = gamma_new / gamma;
          gamma = gamma_new;

          vec_p.axpy(vec_p, vec_z, beta);
          vec_s.axpy(vec_s, vec_Z, beta);

          /// set new defect with our own method, to not use synchronous _set_new_defect method
          status = _update_defect(norm_def_cur->wait());
          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }

      Status _update_defect(const DataType def_cur_norm)
      {
        // increase iteration count
        ++this->_num_iter;

        // save previous defect
        const DataType def_old = this->_def_cur;

        this->_def_cur = def_cur_norm;

        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));

        // plot?
        if(this->_plot)
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " / " << stringify_fp_fix(this->_def_cur / def_old)
            << std::endl;
        }

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_cur))
          return Status::aborted;

        // is diverged?
        if(this->is_diverged())
          return Status::diverged;

        // minimum number of iterations performed?
        if(this->_num_iter < this->_min_iter)
          return Status::progress;

        // is converged?
        if(this->is_converged())
          return Status::success;

        // maximum number of iterations performed?
        if(this->_num_iter >= this->_max_iter)
          return Status::max_iter;

        // check for stagnation?
        if(this->_min_stag_iter > Index(0))
        {
          // did this iteration stagnate?
          if(this->_def_cur >= this->_stag_rate * def_old)
          {
            // increment stagnation count
            if(++this->_num_stag_iter >= this->_min_stag_iter)
              return Status::stagnated;
          }
          else
          {
            // this iteration did not stagnate
            this->_num_stag_iter = Index(0);
          }
        }

        // continue iterating
        return Status::progress;
      }
    }; // class GroppPCG<...>

    /**
     * \brief Creates a new GroppPCG solver object
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
     * A shared pointer to a new GroppPCG object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<GroppPCG<Matrix_, Filter_>> new_gropppcg(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<GroppPCG<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<GroppPCG<Matrix_, Filter_>> new_gropppcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<GroppPCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<GroppPCG<Matrix_, Filter_>> new_gropppcg(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<GroppPCG<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_GROPPPCG_HPP
