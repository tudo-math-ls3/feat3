#pragma once
#ifndef KERNEL_SOLVER_HESSIAN_PRECOND
#define KERNEL_SOLVER_HESSIAN_PRECOND 1
#include <kernel/base_header.hpp>
#include <kernel/solver/nlopt_precond.hpp>
#include <kernel/solver/test_aux/analytic_function_operator.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Preconditioner that applies the inverse of an operator's Hessian
     *
     * \tparam Operator_
     * The operater whose Hessian is used.
     *
     * \tparam Filter_
     * The filter for the preconditioner.
     *
     * \see Solver::NLOptPrecond
     *
     * Note that only specialisations in Operator_ are implemented at the moment.
     *
     * \author Jordi Paul
     */
#ifndef DOXYGEN
    template<typename Operator_, typename Filter_>
    class HessianPrecond;
#else
    template<typename Operator_, typename Filter_>
    class HessianPrecond : public NLOptPrecond<typename Operator_::VectorType, Filter_>
    {};
#endif

    /// \cond internal
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
    template<typename Mem_, typename DT_, typename IT_, typename Function_, typename Filter_>
    class HessianPrecond<AnalyticFunctionOperator<Mem_, DT_, IT_, Function_>, Filter_>
    : public NLOptPrecond<typename AnalyticFunctionOperator<Mem_, DT_, IT_, Function_>::VectorTypeL, Filter_>
    {
      public:
        /// The operator whose Hessian we use for preconditioning
        typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
        /// Our filter type
        typedef Filter_ FilterType;

        /// The input type to apply the Hessian to
        typedef typename OperatorType::VectorTypeL VectorType;
        /// Our base class
        typedef NLOptPrecond<VectorType, Filter_> BaseClass;

        /// Floating point data type
        typedef typename OperatorType::DataType DataType;
        /// Type of the operator's hessian
        typedef typename OperatorType::HessianType HessianType;

      protected:
        /// The operator
        OperatorType& _op;
        /// The filter
        const FilterType& _filter;
        /// We save the operator's Hessian to this
        HessianType _hessian;
        /// And the inverse to this
        HessianType _inv_hessian;

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
          _filter(filter_),
          _hessian(DT_(0)),
          _inv_hessian(DT_(0))
          {
          }

        /**
         * \brief Empty virtual destructor
         */
        virtual ~HessianPrecond()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "HessianPrecond";
        }

        /// \copydoc BaseClass::prepare()
        virtual void prepare(const VectorType& DOXY(vec_state), FilterType& DOXY(filter)) override
        {
          _op.compute_hess(_hessian);
          _inv_hessian.set_inverse(_hessian);
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), _op.name()));

          vec_cor(0, _inv_hessian*vec_def(0));
          this->_filter.filter_cor(vec_cor);

          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 1));

          return Status::success;
        }

        /// \copydoc BaseClass::print()
        virtual void print() const override
        {
          Dist::Comm comm(Dist::Comm::world());
          comm.print(name());
        }
    };
    /// \endcond

    /**
     * \brief Creates a new HessianPrecond
     *
     * \tparam Operator_
     * The operator's type.
     *
     * \tparam Filter_
     * The type of the filte for the preconditioner.
     *
     * \param[in] op_
     * The (nonlinear) operator.
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
     * \brief Wrapper class for applying an the inverse of an operator's approximate Hessian as a preconditioner
     *
     * \tparam Operator_
     * The operator's type
     *
     * \tparam Filter_
     * The type of the filter to apply
     *
     * \see Solver::NLOptPrecond
     *
     * Note that only specialisations in Operator_ are implemented at the moment.
     *
     *
     * \author Jordi Paul
     *
     */
#ifndef DOXYGEN
    template<typename Operator_, typename Filter_>
    class ApproximateHessianPrecond;
#else
    template<typename Operator_, typename Filter_>
    class ApproximateHessianPrecond : public NLOptPrecond<typename OperatorType_::VectorType, Filter>
    {};
#endif

    /// \cond internal
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
    template<typename Mem_, typename DT_, typename IT_, typename Function_, typename Filter_>
    class ApproximateHessianPrecond<AnalyticFunctionOperator<Mem_, DT_, IT_, Function_>, Filter_>
    : public NLOptPrecond<typename AnalyticFunctionOperator<Mem_, DT_, IT_, Function_>::VectorTypeL, Filter_>
    {
      public:
        /// The operator whose Hessian we use for preconditioning
        typedef AnalyticFunctionOperator<Mem_, DT_, IT_, Function_> OperatorType;
        /// Our filter type
        typedef Filter_ FilterType;

        /// The input type to apply the Hessian to
        typedef typename OperatorType::VectorTypeL VectorType;
        /// Our base class
        typedef NLOptPrecond<VectorType, Filter_> BaseClass;

        /// Floating point data type
        typedef typename OperatorType::DataType DataType;
        /// Type of the operator's hessian
        typedef typename OperatorType::HessianType HessianType;

      protected:
        /// The operator
        OperatorType& _op;
        /// The filter
        const FilterType& _filter;
        /// We save the operator's Hessian to this
        HessianType _hessian;
        /// And the inverse to this
        HessianType _inv_hessian;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] op_
         * The operator whose approximate Hessian is to be used.
         *
         * \param[in] filter_
         * The filter to be used.
         *
         */
        explicit ApproximateHessianPrecond(OperatorType& op_, const FilterType& filter_) :
          _op(op_),
          _filter(filter_),
          _hessian(DT_(0)),
          _inv_hessian(DT_(0))
          {
          }

        /**
         * \brief Empty virtual destructor
         *
         */
        virtual ~ApproximateHessianPrecond()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "ApproximateHessianPrecond";
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


        /// \copydoc BaseClass::prepare()
        virtual void prepare(const VectorType& DOXY(vec_state), FilterType& DOXY(filter)) override
        {
          _op.compute_approx_hess(_hessian);
          _inv_hessian.set_inverse(_hessian);
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), _op.name()));

          vec_cor(0, _inv_hessian*vec_def(0));
          this->_filter.filter_cor(vec_cor);

          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 1));

          return Status::success;
        }

        /// \copydoc BaseClass::print()
        virtual void print() const override
        {
          Dist::Comm comm(Dist::Comm::world());
          comm.print(name());
        }
    };
    /// \endcond

    /**
     * \brief Creates a new ApproximateHessianPrecond
     *
     * \tparam Operator_
     * The operator's type.
     *
     * \tparam Filter_
     * The type of the filte for the preconditioner.
     *
     * \param[in] op_
     * The (nonlinear) operator.
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
} // namespace FEAT

#endif // KERNEL_SOLVER_HESSIAN_PRECOND
