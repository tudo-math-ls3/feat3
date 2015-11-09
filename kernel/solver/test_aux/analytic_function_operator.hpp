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
     * Because AnalyticFunctions work in Mem::Main, setting Mem_ to anything else will greatly slow things down due to
     * excessive up-/downloading. It's still here to provide a unified class interface for operators.
     *
     */
    template<typename Mem_, typename DT_, typename IT_, typename Function_>
    class AnalyticFunctionOperator
    {
      public:
        static_assert(std::is_same<typename Function_::ImageType, Analytic::Image::Scalar>::value,
        "AnalyticFunctionOperator is implemented for scalar functions only");

        typedef Mem_ MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, 2> VectorTypeR;
        typedef LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, 2> VectorTypeL;
        typedef VectorTypeR GradientType;

        typedef Analytic::EvalTraits<DT_, Function_> FuncEvalTraits;

        typedef typename Function_::template Evaluator<FuncEvalTraits> EvalType;
        typedef typename EvalType::PointType PointType;

      private:
        Function_& _function;
        EvalType _func_eval;

        IT_ _num_func_evals;
        IT_ _num_grad_evals;
        IT_ _num_hess_evals;

        PointType _my_state;
        typename EvalType::GradientType _my_grad;
        typename EvalType::HessianType _my_hess;
        typename EvalType::HessianType _my_hess_inv;

      public:
        explicit AnalyticFunctionOperator(Function_& func) :
          _function(func),
          _func_eval(func),
          _num_func_evals(IT_(0)),
          _num_grad_evals(IT_(0)),
          _num_hess_evals(IT_(0)),
          _my_state(DT_(0)),
          _my_grad(DT_(0)),
          _my_hess(DT_(0)),
          _my_hess_inv(DT_(0))
          {
          }

        void prepare(const VectorTypeR& vec_state)
        {
          _my_state = vec_state(0);
        }

        typename EvalType::ValueType compute_fval()
        {
          ++_num_func_evals;
          typename EvalType::ValueType ret;
          _func_eval.value(ret, _my_state);
          return ret;
        }

        void compute_gradient(GradientType& vec_out)
        {
          ++_num_grad_evals;
          _func_eval.gradient(_my_grad, _my_state);
          vec_out(0, _my_grad);
        }

        void apply_diag_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++_num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          // Only diagonal of the hessian
          for(int i(0); i < _my_grad.n; ++i)
            _my_grad[i] = vec_in(0)[i]*_my_hess[i][i];

          //DataType det_hess(_my_hess.det());
          //if(det_hess < DataType(0))
          //  std::cout << "Warning: det(Hess) = " << stringify_fp_sci(det_hess) << " < 0!" << std::endl;

          vec_out(0, _my_grad);
        }

        void apply_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++_num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          // Apply exact hessian
          _my_grad = _my_hess*vec_in(0);

          vec_out(0, _my_grad);
        }

        void apply_inv_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          _func_eval.hessian(_my_hess, _my_state);

          _my_hess_inv.set_inverse(_my_hess);
          _my_grad = _my_hess_inv*vec_in(0);

          vec_out(0, _my_grad);
        }

        void apply_approx_inv_hess(VectorTypeR& vec_out, const VectorTypeR& vec_in)
        {
          ++_num_hess_evals;
          _func_eval.hessian(_my_hess, _my_state);

          for(int i(0); i < _my_grad.n; ++i)
            _my_grad[i] = vec_in(0)[i]/Math::abs(_my_hess[i][i]);


          vec_out(0, _my_grad);
        }

        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(IT_(1));
        }

        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(IT_(1));
        }
    };
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_TEST_ANALYTIC_FUNCTION_OPERATOR
