#pragma once
#ifndef KERNEL_SOLVER_RGCR_HPP
#define KERNEL_SOLVER_RGCR_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief (Preconditioned) Recycling Generalized Conjugate Residual solver implementation
     *
     * This class implements a simple Recycling GCR solver.
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see RGCR Paper\cite Benner2011
     *
     * \note The solver internally stores two vectors per iteration, thus may consume a large amount of memory.
     */
    template<
      typename Matrix_,
      typename Filter_>
    class RGCR :
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
      VectorType _vec_r, _vec_p_hat, _vec_q_hat;
      /// descent vectors for later recycling
      std::vector<VectorType> p_list, q_list;

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
      explicit RGCR(const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("RGCR", precond),
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
      explicit RGCR(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("RGCR", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

      virtual String name() const override
      {
        return "RGCR";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _vec_r = this->_system_matrix.create_vector_r();
        _vec_p_hat = this->_system_matrix.create_vector_r();
        _vec_q_hat = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        this->_vec_q_hat.clear();
        this->_vec_p_hat.clear();
        this->_vec_r.clear();
        q_list.clear();
        p_list.clear();
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
        Status st(_apply_intern(vec_cor, vec_def));
        this->plot_summary(st);
        p_list.resize(p_list.size() / 4);
        q_list.resize(q_list.size() / 4);
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
        p_list.resize(p_list.size() / 4);
        q_list.resize(q_list.size() / 4);
        return st;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_p_hat(this->_vec_p_hat);
        VectorType& vec_q_hat(this->_vec_q_hat);

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // no stored descent vectors left, calculate new one
          if (this->_num_iter >= p_list.size())
          {
            // apply preconditioner to defect vector
            if(!this->_apply_precond(vec_p_hat, vec_r, filter))
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            //vec_p_hat.copy(vec_r);
            matrix.apply(vec_q_hat, vec_p_hat);
            filter.filter_def(vec_q_hat);
            for (Index j(0) ; j < this->_num_iter ; ++j)
            {
              auto beta = vec_q_hat.dot(q_list.at(j));
              vec_q_hat.axpy(q_list.at(j), vec_q_hat, -beta);
              vec_p_hat.axpy(p_list.at(j), vec_p_hat, -beta);
            }
            const auto gamma = DataType(1) / vec_q_hat.norm2();
            vec_q_hat.scale(vec_q_hat, gamma);
            q_list.push_back(vec_q_hat.clone());
            vec_p_hat.scale(vec_p_hat, gamma);
            p_list.push_back(vec_p_hat.clone());
          }

          auto alpha = vec_r.dot(q_list.at(this->_num_iter));
          vec_sol.axpy(p_list.at(this->_num_iter), vec_sol, alpha);
          vec_r.axpy(q_list.at(this->_num_iter), vec_r, -alpha);
          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_sol);
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
    }; // class RGCR<...>

    /**
     * \brief Creates a new RGCR solver object
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
     * A shared pointer to a new RGCR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(matrix, filter, precond);
    }
#endif

    /**
     * \brief Creates a new RGCR solver object using a PropertyMap
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
     * A shared pointer to a new RGCR object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(section_name, section, matrix, filter, nullptr);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<RGCR<Matrix_, Filter_>> new_rgcr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<RGCR<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_RGCR_HPP
