#pragma once
#ifndef KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR
#define KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR 1
#include <kernel/base_header.hpp>
#include <kernel/analytic/function.hpp>

#include <deque>

namespace FEAST
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

        /// Class that knows what Function_ wants as input parameters for its evaluator etc.
        typedef Analytic::EvalTraits<DT_, Function_> FuncEvalTraits;
        /// Type of Function_'s evaluator
        typedef typename Function_::template Evaluator<FuncEvalTraits> EvalType;
        /// Input type for Function_
        typedef typename EvalType::PointType PointType;

        static constexpr int dim = PointType::n;

        /// Input vector for the operator and its gradient
        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim> VectorTypeR;
        /// Output vector for the operator's gradient
        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim> VectorTypeL;
        /// Output vector for the operator's gradient
        typedef VectorTypeR GradientType;

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
        /// Temporary variable to save the current inverse hessian
        typename EvalType::HessianType _my_hess_inv;

      public:
        /// Counter for number of function evaluations
        IT_ num_func_evals;
        /// Counter for number of gradient evaluations
        IT_ num_grad_evals;
        /// Counter for number of hessian evaluations
        IT_ num_hess_evals;

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
          _my_hess_inv(DT_(0)),
          num_func_evals(IT_(0)),
          num_grad_evals(IT_(0)),
          num_hess_evals(IT_(0))
          {
          }

        /**
         * \brief Prepares the operator for evaluation by setting the current state
         *
         * \param[in] vec_state
         * The current state
         *
         */
        void prepare(const VectorTypeR& vec_state)
        {
          _my_state = vec_state(0);
        }

        /**
         * \brief Evaluates the operator at the current state
         *
         * \returns The function value at the current state
         */
        typename EvalType::ValueType compute_fval()
        {
          ++num_func_evals;
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
        void compute_gradient(GradientType& vec_out)
        {
          ++num_grad_evals;
          _func_eval.gradient(_my_grad, _my_state);
          vec_out(0, _my_grad);
        }

        /**
         * \brief Applies the diagonal of the hessian at the current state onto a given vector
         *
         * \param[out] vec_out
         * diag(Hess(state)) vec_in
         *
         * \param[in] vec_in
         * (Defect/gradient) vector to apply the hessian to
         */
        void apply_diag_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          // Only diagonal of the hessian
          for(int i(0); i < _my_grad.n; ++i)
            _my_grad[i] = vec_in(0)[i]*_my_hess[i][i];

          //DataType det_hess(_my_hess.det());
          //if(det_hess < DataType(0))
          //  std::cout << "Warning: det(Hess) = " << stringify_fp_sci(det_hess) << " < 0!" << std::endl;

          vec_out(0, _my_grad);
        }

        /**
         * \brief Applies the hessian at the current state onto a given vector
         *
         * \param[out] vec_out
         * Hess(state) vec_in
         *
         * \param[in] vec_in
         * (Defect/gradient) vector to apply the hessian to
         */
        void apply_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          // Apply exact hessian
          _my_grad = _my_hess*vec_in(0);

          vec_out(0, _my_grad);
        }

        /**
         * \brief Applies the inverse hessian at the current state onto a given vector
         *
         * \param[out] vec_out
         * inv(Hess(state)) vec_in
         *
         * \param[in] vec_in
         * (Defect/gradient) vector to apply the inverse hessian to
         */
        void apply_inv_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          _my_hess_inv.set_inverse(_my_hess);
          _my_grad = _my_hess_inv*vec_in(0);

          vec_out(0, _my_grad);
        }

        /**
         * \brief Applies the inverse of the diagonal of the hessian at the current state onto a given vector
         *
         * \param[out] vec_out
         * inv(diag(Hess(state))) vec_in
         *
         * \param[in] vec_in
         * (Defect/gradient) vector to apply the inverse of the diagonal of the hessian to
         */
        void apply_approx_inv_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          for(int i(0); i < _my_grad.n; ++i)
            _my_grad[i] = vec_in(0)[i]/Math::abs(_my_hess[i][i]);

          vec_out(0, _my_grad);
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
    };
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR
