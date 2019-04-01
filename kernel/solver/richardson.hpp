// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_RICHARDSON_HPP
#define KERNEL_SOLVER_RICHARDSON_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Richardson solver implementation
     *
     * This class implements a simple preconditioned Richardson iteration solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class Richardson :
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
      /// damping parameter
      DataType _omega;
      /// defect vector
      VectorType _vec_def;
      /// correction vector
      VectorType _vec_cor;

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
       * \param[in] omega
       * The damping parameter for the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit Richardson(const MatrixType& matrix, const FilterType& filter,
        DataType omega = DataType(1), std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("Richardson", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _omega(omega)
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
      explicit Richardson(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("Richardson", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _omega(1)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
        // Check if we have set omega
        auto omega_p = section->query("omega");
        if(omega_p.second)
        {
          set_omega(DataType(std::stod(omega_p.first)));
        }
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~Richardson()
      {
      }

      virtual String name() const override
      {
        return "Richardson";
      }

      /**
       * \brief Sets the damping parameter
       *
       * \param[in] omega
       * The new damping parameter.
       *
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        _omega = omega;
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        // create two temporary vectors
        _vec_def = this->_system_matrix.create_vector_r();
        _vec_cor = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_cor.clear();
        this->_vec_def.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_def.copy(vec_def);

        // clear solution vector
        vec_cor.format();

        // apply
        this->_status = _apply_intern(vec_cor, vec_def);
        this->plot_summary();
        return this->_status;
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_def);

        // apply
        this->_status = _apply_intern(vec_sol, vec_rhs);
        this->plot_summary();
        return this->_status;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        VectorType& vec_def(this->_vec_def);
        VectorType& vec_cor(this->_vec_cor);
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(vec_def, vec_sol);

        pre_iter.destroy();

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // apply preconditioner
          if(!this->_apply_precond(vec_cor, vec_def, filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }
          //filter.filter_cor(vec_cor);

          // update solution vector
          vec_sol.axpy(vec_cor, vec_sol, this->_omega);

          // compute new defect vector
          matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
          filter.filter_def(vec_def);

          // compute new defect norm
          status = this->_set_new_defect(vec_def, vec_sol);
        }

        // return our status
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
        return status;
      }
    }; // class Richardson<...>

    /**
     * \brief Creates a new Richardson solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The damping parameter for the solver.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new Richardson object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1))
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(matrix, filter, omega, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(matrix, filter, omega, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1),
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(matrix, filter, omega, precond);
    }
#endif

    /**
     * \brief Creates a new Richardson solver object using a PropertyMap
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
     * A shared pointer to a new Richardson object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter, std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Richardson<Matrix_, Filter_>> new_richardson(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<Richardson<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_RICHARDSON_HPP
