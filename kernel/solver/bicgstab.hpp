#pragma once
#ifndef KERNEL_SOLVER_BICGSTAB_HPP
#define KERNEL_SOLVER_BICGSTAB_HPP 1

// includes, FEAST
#include <kernel/solver/iterative.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Bi-Conjugate Gradient Stabilized solver implementation
     *
     * This class implements a simple BiCGStab solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Christoph Lohmann
     */
    template<
      typename Matrix_,
      typename Filter_>
    class BiCGStab :
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
      VectorType _vec_r;
      VectorType _vec_r_tilde;
      VectorType _vec_r_tilde_0;
      VectorType _vec_p_tilde;
      VectorType _vec_v;
      VectorType _vec_v_tilde;
      VectorType _vec_s;
      VectorType _vec_s_tilde;
      VectorType _vec_t;
      VectorType _vec_t_tilde;

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
      explicit BiCGStab(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("BiCGStab", precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      virtual String name() const override
      {
        return "BiCGStab";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create all temporary vectors
        // Each *_tilde vector is a correction, all others are defects.
        _vec_r         = this->_system_matrix.create_vector_r();
        _vec_r_tilde   = this->_system_matrix.create_vector_r();
        _vec_r_tilde_0 = this->_system_matrix.create_vector_r();
        _vec_p_tilde   = this->_system_matrix.create_vector_r();
        _vec_v         = this->_system_matrix.create_vector_r();
        _vec_v_tilde   = this->_system_matrix.create_vector_r();
        _vec_s         = this->_system_matrix.create_vector_r();
        _vec_s_tilde   = this->_system_matrix.create_vector_r();
        _vec_t         = this->_system_matrix.create_vector_r();
        _vec_t_tilde   = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_r.clear();
        this->_vec_r_tilde.clear();
        this->_vec_r_tilde_0.clear();
        this->_vec_p_tilde.clear();
        this->_vec_v.clear();
        this->_vec_v_tilde.clear();
        this->_vec_s.clear();
        this->_vec_s_tilde.clear();
        this->_vec_t.clear();
        this->_vec_t_tilde.clear();
        BaseClass::done_symbolic();
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save rhs vector as initial defect
        this->_vec_r.copy(vec_def);
        //this->_system_filter.filter_def(this->_vec_r);

        // format solution vector
        vec_cor.format();
        return _apply_intern(vec_cor, vec_def);
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        VectorType& vec_r        (this->_vec_r);
        VectorType& vec_r_tilde  (this->_vec_r_tilde);
        VectorType& vec_r_tilde_0(this->_vec_r_tilde_0);
        VectorType& vec_p_tilde  (this->_vec_p_tilde);
        VectorType& vec_v        (this->_vec_v);
        VectorType& vec_v_tilde  (this->_vec_v_tilde);
        VectorType& vec_s        (this->_vec_s);
        VectorType& vec_s_tilde  (this->_vec_s_tilde);
        VectorType& vec_t        (this->_vec_t);
        VectorType& vec_t_tilde  (this->_vec_t_tilde);
        const MatrixType& mat_sys(this->_system_matrix);
        const FilterType& fil_sys(this->_system_filter);
        Status status(Status::progress);

        DataType rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;

        while(status == Status::progress)
        {
          if (restarted == false)
          {
            // initial defect is already computed
            status = this->_set_initial_defect(vec_r, vec_sol);
          }
          else
          {
            mat_sys.apply(vec_r, vec_sol, vec_rhs, -DataType(1));
            fil_sys.filter_def(vec_r);
          }

          // apply preconditioner
          if(!this->_apply_precond(vec_r_tilde_0, vec_r, fil_sys))
            return Status::aborted;
          //fil_sys.filter_cor(vec_r_tilde_0);

          vec_r_tilde.copy(vec_r_tilde_0);
          vec_p_tilde.copy(vec_r_tilde_0);

          rho_tilde = vec_r_tilde_0.dot(vec_r_tilde_0);

          // main BiCGStab loop
          while(status == Status::progress)
          {
            TimeStamp at;
            double mpi_start(Statistics::get_time_mpi_execute());

            mat_sys.apply(vec_v, vec_p_tilde);
            fil_sys.filter_def(vec_v);
            // apply preconditioner
            if(!this->_apply_precond(vec_v_tilde, vec_v, fil_sys))
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              double mpi_stop(Statistics::get_time_mpi_execute());
              Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
              return Status::aborted;
            }
            //fil_sys.filter_cor(vec_v_tilde);

            gamma_tilde = vec_v_tilde.dot(vec_r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              double mpi_stop(Statistics::get_time_mpi_execute());
              Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * vec_v_tilde.norm2()) / this->_def_cur < 1e-5)
            {
              restarted = true;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DataType malpha_tilde(-alpha_tilde);
            vec_s.axpy(vec_v, vec_r, malpha_tilde);

            // compute new defect norm
            /*status = this->_set_new_defect(vec_s);
            if (status == Status::success)
            {
            vec_sol.axpy(vec_p_tilde, vec_sol, alpha_tilde);

            //early_exit = 1;
            //std::cout << "Breakpoint 3 (converged)" << std::endl;
            return status;
            }*/
            vec_s_tilde.axpy(vec_v_tilde, vec_r_tilde, malpha_tilde);

            mat_sys.apply(vec_t, vec_s_tilde);
            fil_sys.filter_def(vec_t);

            // apply preconditioner
            if(!this->_apply_precond(vec_t_tilde, vec_t, fil_sys))
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              double mpi_stop(Statistics::get_time_mpi_execute());
              Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
              return Status::aborted;
            }
            //fil_sys.filter_cor(vec_t_tilde);

            gamma_tilde = vec_t_tilde.dot(vec_t_tilde);
            omega_tilde = vec_t_tilde.dot(vec_s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              double mpi_stop(Statistics::get_time_mpi_execute());
              Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            vec_sol.axpy(vec_s_tilde, vec_sol, omega_tilde);
            vec_sol.axpy(vec_p_tilde, vec_sol, alpha_tilde);

            DataType momega_tilde(-omega_tilde);
            vec_r.axpy(vec_t, vec_s, momega_tilde);

            // compute new defect norm
            status = this->_set_new_defect(vec_r, vec_sol);
            if (status == Status::success)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              double mpi_stop(Statistics::get_time_mpi_execute());
              Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              return status;
            }

            vec_r_tilde.axpy(vec_t_tilde, vec_s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = vec_r_tilde.dot(vec_r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            vec_p_tilde.axpy(vec_v_tilde, vec_p_tilde, momega_tilde);
            vec_p_tilde.scale(vec_p_tilde, beta_tilde);
            vec_p_tilde.axpy(vec_p_tilde, vec_r_tilde);

            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            double mpi_stop(Statistics::get_time_mpi_execute());
            Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
          }
        }

        // finished
        return status;
      }
    }; // class BiCGStab<...>

    /**
     * \brief Creates a new BiCGStab solver object
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
     * A shared pointer to a new BiCGStab object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAST_COMPILER_GNU) && (FEAST_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BiCGStab<Matrix_, Filter_>> new_bicgstab(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<BiCGStab<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_BICGSTAB_HPP
