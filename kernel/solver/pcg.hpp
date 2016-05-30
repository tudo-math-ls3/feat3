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
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class PCG :
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
      /// defect vector
      VectorType _vec_def;
      /// descend direction vector
      VectorType _vec_dir;
      /// temporary vector
      VectorType _vec_tmp;

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

      virtual String name() const override
      {
        return "PCG";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create three temporary vectors
        _vec_def = this->_system_matrix.create_vector_r();
        _vec_dir = this->_system_matrix.create_vector_r();
        _vec_tmp = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_tmp.clear();
        this->_vec_dir.clear();
        this->_vec_def.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_def.copy(vec_def);
        //this->_system_filter.filter_def(this->_vec_def);

        // clear solution vector
        vec_cor.format();

        // apply
        return _apply_intern(vec_cor, vec_def);
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_def);

        // apply
        return _apply_intern(vec_sol, vec_rhs);
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        VectorType& vec_def(this->_vec_def);
        VectorType& vec_dir(this->_vec_dir);
        VectorType& vec_tmp(this->_vec_tmp);
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(vec_def, vec_sol);
        if(status != Status::progress)
          return status;

        // apply preconditioner to defect vector
        if(!this->_apply_precond(vec_dir, vec_def, filter))
          return Status::aborted;
        //filter.filter_cor(vec_dir);

        // compute initial gamma
        DataType gamma = vec_def.dot(vec_dir);

        // start iterating
        while(status == Status::progress)
        {
          TimeStamp at;
          double mpi_start(Statistics::get_time_mpi_execute());

          // compute A*d
          matrix.apply(vec_tmp, vec_dir);
          filter.filter_def(vec_tmp);

          // compute alpha
          DataType alpha = gamma / vec_tmp.dot(vec_dir);

          // update solution vector
          vec_sol.axpy(vec_dir, vec_sol, alpha);

          // update defect vector
          vec_def.axpy(vec_tmp, vec_def, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_def, vec_sol);
          if(status != Status::progress)
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            double mpi_stop(Statistics::get_time_mpi_execute());
            Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
            return status;
          }

          // apply preconditioner
          if(!this->_apply_precond(vec_tmp, vec_def, filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            double mpi_stop(Statistics::get_time_mpi_execute());
            Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
            return Status::aborted;
          }
          //filter.filter_cor(vec_tmp);

          // compute new gamma
          DataType gamma2 = gamma;
          gamma = vec_def.dot(vec_tmp);

          // compute beta
          DataType beta = gamma / gamma2;

          // update direction vector
          vec_dir.axpy(vec_dir, vec_tmp, beta);

          TimeStamp bt;
          Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
          double mpi_stop(Statistics::get_time_mpi_execute());
          Statistics::add_solver_mpi_toe(this->_branch, mpi_stop - mpi_start);
        }

        // we should never reach this point...
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
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_PCG_HPP
