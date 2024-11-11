// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Abstract base class for preconditioners for nonlinear optimization
     *
     * \tparam VectorType_
     * The type of vectors the preconditioner acts upon.
     *
     * \tparam FilterType_
     * The type of filter that gets applied.
     *
     * The main difference to SolverBase is that they have to offer a prepare() function that takes a state variable
     * and prepares the preconditioner. An example would be that prepare() assembles the gradient matrix and apply()
     * applies the inverse of the gradient matrix to the input by using a linear solver.
     *
     * \author Jordi Paul
     *
     */
    template<typename VectorType_, typename FilterType_>
    class NLOptPrecond: public SolverBase<VectorType_>
    {
      public:
        /// Floating point data type
        typedef VectorType_ VectorType;
        /// Type of the filter
        typedef FilterType_ FilterType;

        /**
         * \brief Empty default constructor
         */
        explicit NLOptPrecond()
        {
        }

        /**
         * \brief Empty virtual destructor
         */
        virtual ~NLOptPrecond()
        {
        }

        /**
         * \brief Prepares the preconditioner for application
         *
         * \param[in] vec_state
         * The state variable.
         *
         * \param[in,out] filter
         * The filter that in general will be assembled by this routine.
         *
         */
        virtual void prepare(const VectorType& DOXY(vec_state), FilterType& DOXY(filter))=0;

        /**
         *
         * \brief Returns a descriptive String
         *
         * \returns The class name as String.
         */
        virtual String name() const override = 0;

        /**
         * \brief Applies the preconditioner
         *
         * \param[in, out] vec_cor
         * This holds the initial guess and will be overwritten by the solution.
         *
         * \param[in] vec_def
         * The right hand side the preconditioner is applied to.
         *
         * Usually this will apply some kind of solver.
         *
         * \returns The Solver::Status of the application of the preconditioner.
         */
        virtual Status apply(VectorType& DOXY(vec_cor), const VectorType& DOXY(vec_def)) override = 0;

        /**
         * \brief Returns information
         */
        virtual String info() const
        {
          return this->name();
        }

    }; // class NLOptPrecond

    /**
     * \brief Wrapper class around a (potentially nonlinear) operator
     *
     * \tparam NonlinearOperator_
     * The operator whose inverse is applied.
     *
     * This will apply the inverse of a given operator, i.e. as a preconditioner.
     *
     * \author Jordi Paul
     *
     */
    template<typename NonlinearOperator_>
    class NonlinearOperatorPrecondWrapper
    : public Solver::NLOptPrecond
    <
      typename NonlinearOperator_::SystemLevelType::GlobalSystemVectorR,
      typename NonlinearOperator_::SystemLevelType::GlobalSystemFilter
    >
    {
      public:
        /// The vector type the preconditioner can be applied to
        typedef typename NonlinearOperator_::SystemLevelType::GlobalSystemVectorR VectorType;
        /// The filter type
        typedef typename NonlinearOperator_::SystemLevelType::GlobalSystemFilter FilterType;
        /// Our base class
        typedef Solver::NLOptPrecond<VectorType, FilterType> BaseClass;

      private:
        /// The operator that defines the preconditioner
        NonlinearOperator_ _op;

      public:
        /**
         * \brief Variadic template constructor
         */
        template<typename... Args_>
        explicit NonlinearOperatorPrecondWrapper(Args_&&... args):
          _op(std::forward<Args_>(args)...)
          {
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~NonlinearOperatorPrecondWrapper()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "NonlinearOperatorPrecondWrapper<"+_op.name()+">";
        }

        /// \copydoc BaseClass::prepare()
        virtual void prepare(const VectorType& vec_state, FilterType& DOXY(filter)) override
        {
          _op.prepare(vec_state);
        }

        /// \copydoc BaseClass::init_numeric()
        virtual void init_numeric() override
        {
          _op.init_numeric();
        }

        /// \copydoc BaseClass::apply()
        virtual Solver::Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), _op.name()));

          Solver::Status st(_op.apply(vec_cor, vec_def));

          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), st, 1));

          return st;
        }

        /// \copydoc BaseClass::info()
        virtual String info() const override
        {
          return this->name() + "[" + _op.info() + "]";
        }

    }; // class NonlinearOperatorPrecondWrapper

    /**
     * \brief Creates a new NonlinearOperatorPrecondWrapper
     *
     * \tparam NonlinearOperator_
     * The operator to be used.
     *
     * \param[in] args
     * Arguments passed to the NonlinearOperator's constructor.
     *
     * \returns A shared_ptr to a new NonlinearOperatorPrecondWrapper object.
     */
    template<typename NonlinearOperator_, typename... Args_>
    inline std::shared_ptr<NonlinearOperatorPrecondWrapper<NonlinearOperator_>>
    new_nonlinear_operator_precond_wrapper(Args_&&... args)
    {
      return std::make_shared<NonlinearOperatorPrecondWrapper<NonlinearOperator_>>(std::forward<Args_>(args)...);
    }

  } // namespace Solver
}// namespace FEAT
