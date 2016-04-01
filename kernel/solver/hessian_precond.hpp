#pragma once
#ifndef KERNEL_SOLVER_HESSIAN_PRECOND
#define KERNEL_SOLVER_HESSIAN_PRECOND 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Wrapper class for applying an the inverse of an operator's Hessian as a preconditioner
     *
     * \tparam Operator_
     * The operator's type
     *
     * \tparam Filter_
     * The type of the filter to apply
     *
     * \author Jordi Paul
     *
     */
    template<typename Operator_, typename Filter_>
    class HessianPrecond : public SolverBase<typename Operator_::VectorTypeL>
    {
      public:
        /// The operator whose Hessian we use for preconditioning
        typedef Operator_ OperatorType;
        /// Our filter type
        typedef Filter_ FilterType;
        /// The input type to apply the Hessian to
        typedef typename OperatorType::VectorTypeL VectorType;
        /// Floating point data type
        typedef typename OperatorType::DataType DataType;

      protected:
        /// The operator
        OperatorType& _op;
        /// The filter
        const FilterType& _filter;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] op_
         * The operator whose Hessian is to be used.
         *
         * \param[in] filter_
         * The filter to be used.
         *
         */
        explicit HessianPrecond(OperatorType& op_, const FilterType& filter_) :
          _op(op_),
          _filter(filter_)
          {
          }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return "Hessian";
        }

        //virtual void init_symbolic() override
        //{
        //  _inv_diag = _op.create_vector_r();
        //}

        //virtual void done_symbolic() override
        //{
        //  _inv_diag.clear();
        //}

        //virtual void init_numeric() override
        //{
        //  // extract matrix diagonal
        //  _op.extract_diag(_inv_diag);

        //  // invert diagonal elements
        //  _inv_diag.component_invert(_inv_diag, _omega);
        //}

        /**
         * \brief Applied the preconditioner
         *
         * \param[out] vec_cor
         * The returned preconditioned input
         *
         * \param[in] vec_def
         * The vector to apply the preconditioner to
         *
         * \returns Status::success if the inverse hessian was applied successfully.
         */
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          Status st(_op.apply_inv_hess(vec_cor, vec_def));
          this->_filter.filter_cor(vec_cor);

          return st;
        }
    };

    /**
     * \brief Creates a new HessianPrecond
     *
     * \param[in] op_
     * The (nonlinear) operator
     *
     * \param[in] filter_
     * The system filter.
     *
     * \returns
     * A shared pointer to a new HessianPrecond object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<HessianPrecond<Operator_, Filter_>> new_hessian_precond(
      Operator_& op_, const Filter_& filter_)
      {
        return std::make_shared<HessianPrecond<Operator_, Filter_>>(op_, filter_);
      }

    /**
     * \brief Wrapper class for applying an approximation to the inverse of an operator's Hessian as a preconditioner
     *
     * \tparam Operator_
     * The operator's type
     *
     * \tparam Filter_
     * The type of the filter to apply
     *
     * \author Jordi Paul
     *
     * This nearly does the same as the HessianPrecond to the letter, but inheriting from it would not help much.
     *
     */
    template<typename Operator_, typename Filter_>
    class ApproximateHessianPrecond : public SolverBase<typename Operator_::VectorTypeL>
    {
      public:
        /// The operator whose Hessian we use for preconditioning
        typedef Operator_ OperatorType;
        /// Our filter type
        typedef Filter_ FilterType;
        /// The input type to apply the Hessian to
        typedef typename OperatorType::VectorTypeL VectorType;
        /// Floating point data type
        typedef typename OperatorType::DataType DataType;

      protected:
        /// The operator
        OperatorType& _op;
        /// The filter
        const FilterType& _filter;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] op_
         * The operator whose Hessian is to be used.
         *
         * \param[in] filter_
         * The filter to be used.
         *
         */
        explicit ApproximateHessianPrecond(OperatorType& op_, const FilterType& filter_) :
          _op(op_),
          _filter(filter_)
          {
          }

        /// Returns the name of the solver.
        virtual String name() const override
        {
          return "ApproximateHessian";
        }

        //virtual void init_symbolic() override
        //{
        //  _inv_diag = _op.create_vector_r();
        //}

        //virtual void done_symbolic() override
        //{
        //  _inv_diag.clear();
        //}

        //virtual void init_numeric() override
        //{
        //  // extract matrix diagonal
        //  _op.extract_diag(_inv_diag);

        //  // invert diagonal elements
        //  _inv_diag.component_invert(_inv_diag, _omega);
        //}

        /**
         * \brief Applied the preconditioner
         *
         * \param[out] vec_cor
         * The returned preconditioned input
         *
         * \param[in] vec_def
         * The vector to apply the preconditioner to
         *
         * \returns Status::success if the approximate inverse hessian was applied successfully.
         */
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          Status st(_op.apply_approx_inv_hess(vec_cor, vec_def));
          this->_filter.filter_cor(vec_cor);

          return st;
        }
    };

    /**
     * \brief Creates a new ApproximateHessianPrecond
     *
     * \param[in] op_
     * The (nonlinear) operator
     *
     * \param[in] filter_
     * The system filter.
     *
     * \returns
     * A shared pointer to a new HessianPrecond object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<ApproximateHessianPrecond<Operator_, Filter_>> new_approximate_hessian_precond(
      Operator_& op_, const Filter_& filter_)
      {
        return std::make_shared<ApproximateHessianPrecond<Operator_, Filter_>>(op_, filter_);
      }

  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_HESSIAN_PRECOND
