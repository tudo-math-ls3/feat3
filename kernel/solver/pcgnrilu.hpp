#pragma once
#ifndef KERNEL_SOLVER_PCGNRILU_HPP
#define KERNEL_SOLVER_PCGNRILU_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/ilu_precond.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief ILU(p)-preconditioned Conjugate-Gradient on Normal Equations solver implementation
     *
     * Let \f$M = LU\f$ denote the ILU(p)-factorisation of a (possibly unsymmetric and indefinite)
     * system matrix \f$A\f$, then this class implements a Conjugate-Gradient solver applied onto
     * the preconditioned normal equation system:
     *
     * \f[ \underbrace{(U^{-\top}A^\top L^{-\top} L^{-1}AU^{-1})}_{\widetilde{A}}
     * ~\cdot~\underbrace{(U x)}_{\widetilde{x}}~=~
     * \underbrace{(U^{-\top}A^\top L^{-\top} L^{-1} b)}_{\widetilde{b}}\f]
     *
     * \attention
     * This solver implementation is experimental. Do not use it unless you know exactly what you are doing!
     *
     * \note
     * This class implements a special case of the more generic PCGNR solver with the preconditioners
     * \f$M_L = LL^\top\f$ and \f$M_R = U^\top U\f$.
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
    class PCGNRILU :
      public IterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef typename MatrixType::IndexType IndexType;
      typedef IterativeSolver<VectorType> BaseClass;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// the transposed system matrix
      MatrixType _transp_matrix;
      /// temporary vectors
      VectorType _vec_r, _vec_p, _vec_q, _vec_s, _vec_t;
      /// ilu fill-in
      int _ilu_p;
      /// ilu factorisation data
      Intern::ILUCore<DataType, IndexType> _ilu;

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
       * \param[in] ilu_p
       * Maximum allowed fill-in for ILU preconditioner. Set to -1 to disable preconditioning.
       */
      explicit PCGNRILU(const MatrixType& matrix, const FilterType& filter, int ilu_p = -1) :
        BaseClass("PCGNRILU"),
        _system_matrix(matrix),
        _system_filter(filter),
        _ilu_p(ilu_p)
      {
      }

      virtual String name() const override
      {
        return "PCGNRILU";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        _vec_p = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_s = this->_system_matrix.create_vector_r();
        _vec_t = this->_system_matrix.create_vector_r();

        if(_ilu_p >= 0)
        {
          _ilu.set_struct(this->_system_matrix);
          _ilu.factorise_symbolic(_ilu_p);
          _ilu.alloc_data();
        }
      }

      virtual void done_symbolic() override
      {
        _ilu.clear();

        this->_vec_t.clear();
        this->_vec_s.clear();
        this->_vec_r.clear();
        this->_vec_q.clear();
        this->_vec_p.clear();

        BaseClass::done_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        _transp_matrix.transpose(_system_matrix);
        if(_ilu_p >= 0)
        {
          _ilu.copy_data(_system_matrix);
          _ilu.factorise_numeric_il_du();
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_r.copy(vec_def);

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
      void _precond_l(VectorType& vec_c, const VectorType& vec_q)
      {
        vec_c.copy(vec_q);
        if(_ilu_p >= 0)
        {
          DataType* x = vec_c.elements();
          _ilu.solve_il(x, x);
          _ilu.solve_ilt(x, x);
        }
      }

      void _precond_r(VectorType& vec_c, const VectorType& vec_q)
      {
        vec_c.copy(vec_q);
        if(_ilu_p >= 0)
        {
          DataType* x = vec_c.elements();
          _ilu.solve_dut(x, x);
          _ilu.solve_du(x, x);
        }
      }

      virtual Status _apply_intern(VectorType& vec_x, const VectorType& DOXY(vec_b))
      {
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

        // compute initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_x);
        if(status != Status::progress)
          return status;

        // apply left preconditioner to defect vector
        // p[0] := (L*L^T)^{-1} r[0]
        this->_precond_l(vec_p, vec_r);
        filter.filter_cor(vec_p);

        // apply transposed matrix
        // s[0] := A^T * p[0]
        transp.apply(vec_s, vec_p);
        filter.filter_def(vec_s);

        // apply right preconditioner
        // q[0] := (R^T*R)^{-1} * s[0]
        this->_precond_r(vec_q, vec_s);
        filter.filter_cor(vec_q);

        // compute initial gamma:
        // gamma[0] := < s[0], q[0] >
        DataType gamma = vec_s.dot(vec_q);

        // start iterating
        while(status == Status::progress)
        {
          TimeStamp at;
          double mpi_start(Statistics::get_time_mpi_execute());

          // y[k] := A * q[k]
          matrix.apply(vec_y, vec_q);
          filter.filter_def(vec_y);

          // z[k] := (L*L^T)^{-1} * y[k]
          this->_precond_l(vec_z, vec_y);
          filter.filter_cor(vec_z);

          // compute alpha
          // alpha[k] := gamma[k] / < y[k], z[k] >
          DataType alpha = gamma / vec_y.dot(vec_z);

          // update solution vector
          // x[k+1] := x[k] + alpha[k] * q[k]
          vec_x.axpy(vec_q, vec_x, alpha);

          // update defect vector
          // r[k+1] := r[k] - alpha[k] * y[k]
          vec_r.axpy(vec_y, vec_r, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_x);
          if(status != Status::progress)
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            double mpi_stop(Statistics::get_time_mpi_execute());
            Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
            return status;
          }

          // p[k+1] := p[k] - alpha[k] * z[k]
          vec_p.axpy(vec_z, vec_p, -alpha);

          // s[k+1] := A^T * p[k+1]
          transp.apply(vec_s, vec_p);
          filter.filter_def(vec_s);

          // t[k+1] := (R^T*R)^{-1} * s[k+1]
          this->_precond_r(vec_t, vec_s);
          filter.filter_cor(vec_t);

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

          TimeStamp bt;
          Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
          double mpi_stop(Statistics::get_time_mpi_execute());
          Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
        }

        // we should never reach this point...
        return Status::undefined;
      }
    }; // class PCGNRILU<...>

    /**
     * \brief Creates a new PCGNRILU solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] ilu_p
     * Maximum allowed fill-in for ILU preconditioner. Set to -1 to disable preconditioning.
     *
     * \returns
     * A shared pointer to a new PCGNRILU object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<PCGNRILU<Matrix_, Filter_>> new_pcgnrilu(
      const Matrix_& matrix, const Filter_& filter, int ilu_p = -1)
    {
      return std::make_shared<PCGNRILU<Matrix_, Filter_>>(matrix, filter, ilu_p);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PCGNRILU_HPP
