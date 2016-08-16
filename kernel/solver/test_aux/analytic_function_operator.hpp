#pragma once
#ifndef KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR
#define KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR 1
#include <kernel/base_header.hpp>
#include <kernel/analytic/function.hpp>

#include <deque>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Wrapper class defining an operator based on a scalar AnalyticFunction
     *
     * \tparam Mem_
     * Memory type the optimiser has to work in
     *
     * \tparam DT_
     * Floating point precision
     *
     * \tparam IT_
     * Index type
     *
     * \tparam
     * The scalar AnalyticFunction
     *
     * Since all solvers work on LAFEM containers, this is hardwired to expect LAFEM containers and just use the
     * first entry, i.e. a scalar function with 3 input parameters results in a LAFEM::DenseVectorBlocked with one
     * entry of BlockSize 3.
     *
     * Because AnalyticFunctions work in Mem::Main, setting Mem_ to anything else will greatly slow things down due
     * to excessive up-/downloading. It's still here to provide a unified class interface for operators.
     *
     */
    template<typename Mem_, typename DT_, typename IT_, typename Function_>
    class AnalyticFunctionOperator
    {
      public:
        /// Check if Function_ is scalar
        static_assert(std::is_same<typename Function_::ImageType, Analytic::Image::Scalar>::value,
        "AnalyticFunctionOperator is implemented for scalar functions only");

        /// Memory architecture of the vectors we work in
        typedef Mem_ MemType;
        /// Floating point type
        typedef DT_ DataType;
        /// Integer type
        typedef IT_ IndexType;

        /// The type of the underlying analytic function
        typedef Function_ FunctionType;

        /// Can the operator compute its full Hessian?
        static constexpr bool can_hess = FunctionType::can_hess;
        /// Can the operator compute the diagonal of its Hessian?
        static constexpr bool can_diag_hess = can_hess;

        /// Class that knows what Function_ wants as input parameters for its evaluator etc.
        typedef Analytic::EvalTraits<DT_, Function_> FuncEvalTraits;
        /// Type of Function_'s evaluator
        typedef typename Function_::template Evaluator<FuncEvalTraits> EvalType;
        /// Input type for Function_
        typedef typename EvalType::PointType PointType;

        /// Dimension each variable is from
        static constexpr int dim = PointType::n;

        /// Input vector for the operator and its gradient
        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim> VectorTypeR;
        /// Output vector for the operator's gradient
        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim> VectorTypeL;
        /// Output vector for the operator's gradient
        typedef VectorTypeR GradientType;
        /// Output matrix for the operator's hessian
        typedef typename EvalType::HessianType HessianType;

      private:
        /// The Analytic::Function
        Function_& _function;
        /// _function's Evaluator
        EvalType _func_eval;

        /// The current point, the function and its gradient/hessian are evaluated here
        PointType _my_state;
        /// Temporary variable to save the current gradient or other vectors
        typename EvalType::GradientType _my_grad;
        /// Temporary variable to save the current hessian or other matrices
        typename EvalType::HessianType _my_hess;

        /// Counter for number of function evaluations
        Index _num_func_evals;
        /// Counter for number of gradient evaluations
        Index _num_grad_evals;
        /// Counter for number of hessian evaluations
        Index _num_hess_evals;

      public:
        /**
         * \brief The one and only sensible constructor
         *
         * \param[in] func
         * The AnalyticFunction defining the operator's action
         *
         */
        explicit AnalyticFunctionOperator(Function_& func) :
          _function(func),
          _func_eval(func),
          _my_state(DT_(0)),
          _my_grad(DT_(0)),
          _my_hess(DT_(0)),
          _num_func_evals(Index(0)),
          _num_grad_evals(Index(0)),
          _num_hess_evals(Index(0))
          {
          }

        static String name()
        {
          return "AnalyticFunctionOperator";
        }

        /**
         * \brief Creates an empty L-vector of appropriate size
         *
         * \returns Empty L-vector of length 1
         */
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(IT_(1));
        }

        /**
         * \brief Creates an empty R-vector of appropriate size
         *
         * \returns Empty R-vector of length 1
         */
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(IT_(1));
        }

        /**
         * \brief The number of input variables for the operator and its gradient
         *
         * For compatibility with linear operators represented by matrices, this is called columns.
         *
         * \returns The number of input variables the operator needs.
         */
        Index columns()
        {
          return Index(dim);
        }

        /**
         * \brief The number of output variables of the operator's gradient
         *
         * For compatibility with linear operators represented by matrices, this is called rows.
         *
         * \returns The number of out variables the operator's gradient gives.
         */
        Index rows()
        {
          return Index(dim);
        }

        /**
         * \returns The number of times the functional value was computed.
         */
        Index get_num_func_evals() const
        {
          return _num_func_evals;
        }

        /**
         * \returns The number of times the gradient was computed.
         */
        Index get_num_grad_evals() const
        {
          return _num_grad_evals;
        }

        /**
         * \returns The number of times the hessian was computed.
         */
        Index get_num_hess_evals() const
        {
          return _num_hess_evals;
        }

        /**
         * \brief Resets all evaluation counts
         */
        void reset_num_evals()
        {
          _num_func_evals = Index(0);
          _num_grad_evals = Index(0);
          _num_hess_evals = Index(0);
        }

        /**
         * \brief Prepares the operator for evaluation by setting the current state
         *
         * \param[in] vec_state
         * The current state.
         *
         * \param[in] filter
         * The system filter. Kept for a unified interface because other operators need to re-assemble it.
         *
         */
        template<typename FilterType_>
        void prepare(const VectorTypeR& vec_state, FilterType_ DOXY(filter))
        {
          _my_state = vec_state(0);
        }

        /**
         * \brief Evaluates the operator at the current state
         *
         * \returns The function value at the current state
         */
        typename EvalType::ValueType compute_func()
        {
          ++_num_func_evals;
          typename EvalType::ValueType ret;
          _func_eval.value(ret, _my_state);
          return ret;
        }

        /**
         * \brief Evaluates the gradient of the operator at the current state
         *
         * \param[out] vec_out
         * The vector containing the gradient
         */
        void compute_grad(GradientType& vec_out)
        {
          ++_num_grad_evals;
          _func_eval.gradient(_my_grad, _my_state);
          vec_out(0, _my_grad);
        }

        /**
         * \brief Applies the Hessian of the operator at the current state to a vector
         *
         * \param[out] vec_out
         * Hessian*vec_in
         *
         * \param[in] vec_in
         * Vector the Hessian gets applied to.
         */
        void apply_hess(VectorTypeL& vec_out, const VectorTypeR& vec_in)
        {
          ++_num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);
          vec_out(0, _my_hess*vec_in(0));
        }

        /**
         * \brief Evaluates the Hessian of the operator at the current state
         *
         * \param[out] mat_out
         * The matrix containing the hessiant
         */
        void compute_hess(HessianType& mat_out)
        {
          ++_num_hess_evals;
          _func_eval.hessian(mat_out, _my_state);
        }

        /**
         * \brief Evaluates the Hessian of the operator at the current state
         *
         * \param[out] mat_out
         * The matrix containing the hessiant
         */
        void compute_approx_hess(HessianType& mat_out)
        {
          ++_num_hess_evals;
          _func_eval.hessian(mat_out, _my_state);
          for(int i(1); i < mat_out.m; ++i)
          {
            for(int j(0); j < i; ++j)
            {
              mat_out(i,j) = DT_(0);
              mat_out(j,i) = DT_(0);
            }
          }
        }
    };
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR
