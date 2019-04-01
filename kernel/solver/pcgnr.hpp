// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_PCGNR_HPP
#define KERNEL_SOLVER_PCGNR_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Conjugate-Gradient on Normal Equations solver implementation
     *
     * Let \f$M_L = LL^\top\f$ and \f$M_R = R^\top R\f$ denote two <b>symmetric and positive definite</b>
     * preconditioners for the (possibly unsymmetric and indefinite) system matrix \f$A\f$, then this class
     * implements a Conjugate-Gradient solver applied onto the preconditioned normal equation system:
     *
     * \f[ \underbrace{(R^{-\top}A^\top L^{-\top} L^{-1}AR^{-1})}_{\widetilde{A}}
     * ~\cdot~\underbrace{(R x)}_{\widetilde{x}}~=~
     * \underbrace{(R^{-\top}A^\top L^{-\top} L^{-1} b)}_{\widetilde{b}}\f]
     *
     * \attention
     * This solver implementation is experimental. Do not use it unless you know exactly what you are doing!
     *
     * \note
     * In contrast to most other iterative solver variants, this class does not derive from
     * PreconditionedIterativeSolver, as there are two possibly different preconditioners in here.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PCGNR :
      public IterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef IterativeSolver<VectorType> BaseClass;

      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// left preconditioner
      std::shared_ptr<PrecondType> _precond_l;
      /// right preconditioner
      std::shared_ptr<PrecondType> _precond_r;
      /// the transposed system matrix
      MatrixType _transp_matrix;
      /// temporary vectors
      VectorType _vec_r, _vec_p, _vec_q, _vec_s, _vec_t;

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
       * \param[in] precond_l
       * A pointer to the left preconditioner M_L. May be \c nullptr.
       *
       * \param[in] precond_r
       * A pointer to the right preconditioner M_R. May be \c nullptr.
       *
       * \note
       * \p precond_l and \p precond_r may point to the same object.
       */
      explicit PCGNR(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond_l = nullptr,
        std::shared_ptr<PrecondType> precond_r = nullptr) :
        BaseClass("PCGNR"),
        _system_matrix(matrix),
        _system_filter(filter),
        _precond_l(precond_l),
        _precond_r(precond_r)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] precond_l
       * A pointer to the left preconditioner M_L. May be \c nullptr.
       *
       * \param[in] precond_r
       * A pointer to the right preconditioner M_R. May be \c nullptr.
       *
       * \note
       * \p precond_l and \p precond_r may point to the same object.
       */
      explicit PCGNR(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond_l = nullptr,
        std::shared_ptr<PrecondType> precond_r = nullptr) :
        BaseClass("PCGNR", section_name, section),
        _system_matrix(matrix),
        _system_filter(filter),
        _precond_l(precond_l),
        _precond_r(precond_r)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
      }

      virtual String name() const override
      {
        return "PCGNR";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // initialise preconditioners
        if(_precond_l)
          _precond_l->init_symbolic();
        if((_precond_r) && (_precond_r != _precond_l))
          _precond_r->init_symbolic();

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

        // release preconditioners
        if((_precond_r) && (_precond_r != _precond_l))
          _precond_r->done_symbolic();
        if(_precond_l)
          _precond_l->done_symbolic();

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // initialise preconditioners
        if(_precond_l)
          _precond_l->init_numeric();
        if((_precond_r) && (_precond_r != _precond_l))
          _precond_r->init_numeric();

        // transpose matrix
        this->_transp_matrix = this->_system_matrix.transpose();
      }

      virtual void done_numeric() override
      {
        // release transposed matrix
        this->_transp_matrix.clear();

        // release preconditioners
        if((_precond_r) && (_precond_r != _precond_l))
          _precond_r->done_numeric();
        if(_precond_l)
          _precond_l->done_numeric();

        BaseClass::done_numeric();
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
      bool _apply_precond_l(VectorType& vec_cor, const VectorType& vec_def)
      {
        if(_precond_l)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_precond_l->name()));
          return status_success(_precond_l->apply(vec_cor, vec_def));
        }
        vec_cor.copy(vec_def);
        this->_system_filter.filter_cor(vec_cor);
        return true;
      }

      bool _apply_precond_r(VectorType& vec_cor, const VectorType& vec_def)
      {
        if(_precond_r)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_precond_r->name()));
          return status_success(_precond_r->apply(vec_cor, vec_def));
        }
        vec_cor.copy(vec_def);
        this->_system_filter.filter_cor(vec_cor);
        return true;
      }

      virtual Status _apply_intern(VectorType& vec_sol)
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const MatrixType& transp(this->_transp_matrix);
        const FilterType& filter(this->_system_filter);

        VectorType& vec_p(this->_vec_p);
        VectorType& vec_q(this->_vec_q);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_s(this->_vec_s);
        VectorType& vec_t(this->_vec_t);
        VectorType& vec_y(this->_vec_s); // same as s
        VectorType& vec_z(this->_vec_t); // same as t

        // set initial defect
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        // apply left preconditioner to defect vector
        // p[0] := M_L^{-1} * r[0]
        if(!this->_apply_precond_l(vec_p, vec_r))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        // apply transposed matrix
        // s[0] := A^T * p[0]
        transp.apply(vec_s, vec_p);

        // apply right preconditioner
        // q[0] := M_R^{-1} * s[0]
        if(!this->_apply_precond_r(vec_q, vec_s))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        // compute initial gamma
        // gamma[0] := < s[0], q[0] >
        DataType gamma = vec_s.dot(vec_q);

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // y[k] := A * q[k]
          matrix.apply(vec_y, vec_q);
          filter.filter_def(vec_y);

          // z[k] := M_L^{-1} * y[k]
          if(!this->_apply_precond_l(vec_z, vec_y))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // compute alpha
          // alpha[k] := gamma[k] / < y[k], z[k] >
          DataType alpha = gamma / vec_y.dot(vec_z);

          // update solution vector
          // x[k+1] := x[k] + alpha[k] * q[k]
          vec_sol.axpy(vec_q, vec_sol, alpha);

          // update defect vector
          // r[k+1] := r[k] - alpha[k] * y[k]
          vec_r.axpy(vec_y, vec_r, -alpha);

          // compute defect norm and check for convergence
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // p[k+1] := p[k] - alpha[k] * z[k]
          vec_p.axpy(vec_z, vec_p, -alpha);

          // s[k+1] := A^T * p[k+1]
          transp.apply(vec_s, vec_p);
          filter.filter_def(vec_s);

          // t[k+1] := M_R^{-1} * s[k+1]
          if(!this->_apply_precond_r(vec_t, vec_s))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // compute new gamma
          // gamma[k+1] := < s[k+1], t[k+1] >
          DataType gamma2 = gamma;
          gamma = vec_s.dot(vec_t);

          // compute beta
          // beta[k] := gamma[k+1] / gamma[k]
          DataType beta = gamma / gamma2;

          // update direction vector
          // q[k+1] := t[k+1] + beta * q[k]
          vec_q.axpy(vec_q, vec_t, beta);
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }
    }; // class PCGNR<...>

    /**
     * \brief Creates a new PCGNR solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond_l
     * The left preconditioner. May be \c nullptr.
     *
     * \param[in] precond_r
     * The right preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new PCGNR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(matrix, filter);
    }

    template<typename Matrix_, typename Filter_, typename PrecondL_, typename PrecondR_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<PrecondL_> precond_l,
      std::shared_ptr<PrecondR_> precond_r)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(matrix, filter, precond_l, precond_r);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond_l = nullptr,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond_r = nullptr)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(matrix, filter, precond_l, precond_r);
    }
#endif

    /**
     * \brief Creates a new PCGNR solver object using a PropertyMap
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
     * \param[in] precond_l
     * The left preconditioner. May be \c nullptr.
     *
     * \param[in] precond_r
     * The right preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new PCGNR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    template<typename Matrix_, typename Filter_, typename PrecondL_, typename PrecondR_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<PrecondL_> precond_l,
      std::shared_ptr<PrecondR_> precond_r)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond_l, precond_r);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCGNR<Matrix_, Filter_>> new_pcgnr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond_l = nullptr,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond_r = nullptr)
    {
      return std::make_shared<PCGNR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond_l, precond_r);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PCGNR_HPP
