#pragma once
#ifndef KERNEL_ANALYTIC_WRAPPERS_HPP
#define KERNEL_ANALYTIC_WRAPPERS_HPP 1

// includes, FEAT
#include <kernel/analytic/function.hpp>

namespace FEAT
{
  namespace Analytic
  {
    /**
     * \brief Analytic Function Gradient wrapper
     *
     * This class template represents the gradient of another scalar function.
     *
     * \tparam Function_
     * The function whose gradient is to be wrapped.
     *
     * \author Peter Zajac
     */
    template<typename Function_>
    class Gradient :
      public Analytic::Function
    {
    public:
      // ensure that the input function is scalar
      static_assert(Function_::ImageType::is_scalar, "only scalar functions are supported");

      /// our domain dimension is the same as the input function's
      static constexpr int domain_dim = Function_::domain_dim;

      /// this is a vector-valued function
      typedef Image::Vector<domain_dim> ImageType;

      static constexpr bool can_value = Function_::can_grad;
      static constexpr bool can_grad = Function_::can_hess;
      static constexpr bool can_hess = false;

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;
        typedef typename Traits_::HessianType HessianType;

      private:
        /// our original function evaluator
        typename Function_::template Evaluator<EvalTraits<DataType, Function_>> _func_eval;

      public:
        explicit Evaluator(const Gradient& function) :
          _func_eval(function._function)
        {
        }

        void value(ValueType& val, const PointType& point)
        {
          _func_eval.gradient(val, point);
        }

        void gradient(GradientType& grad, const PointType& point)
        {
          _func_eval.hessian(grad, point);
        }
      };

    private:
      const Function_& _function;

    public:
      explicit Gradient(const Function_& function) :
        _function(function)
      {
      }
    };

    /**
     * \brief Analytic Function Divergence wrapper
     *
     * This class template represents the divergence of another vector field.
     *
     * \tparam Function_
     * The function whose divergence is to be wrapped.
     *
     * \author Peter Zajac
     */
    template<typename Function_>
    class Divergence :
      public Analytic::Function
    {
    public:
      // ensure that the input function is scalar
      static_assert(Function_::ImageType::is_vector, "only vector fields are supported");

      /// our domain dimension is the same as the input function's
      static constexpr int domain_dim = Function_::domain_dim;

      /// this is a vector-valued function
      typedef Image::Scalar ImageType;

      static constexpr bool can_value = Function_::can_grad;
      static constexpr bool can_grad = false; //Function_::can_hess;
      static constexpr bool can_hess = false;

      template<typename Traits_>
      class Evaluator :
        public Function::Evaluator<Traits_>
      {
      public:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;
        typedef typename Traits_::HessianType HessianType;

      private:
        typedef EvalTraits<DataType, Function_> FuncEvalTraits;
        /// our original function evaluator
        typename Function_::template Evaluator<FuncEvalTraits> _func_eval;

      public:
        explicit Evaluator(const Divergence& function) :
          _func_eval(function._function)
        {
        }

        void value(ValueType& val, const PointType& point)
        {
          typename FuncEvalTraits::GradientType grad;
          _func_eval.gradient(grad, point);
          val = grad.trace();
        }

        /// \todo implement gradient
        /*void gradient(GradientType& grad, const PointType& point)
        {
          typename FuncEvalTraits::HessianType hess;
          _func_eval.hessian(hess, point);
          for(int i(0); i < domain_dim; ++i)
          {
            grad[i] = DataType(0);
            for(int j()
          }
        }*/
      };

    private:
      const Function_& _function;

    public:
      explicit Divergence(const Function_& function) :
        _function(function)
      {
      }
    };

    /**
     * \brief Analytic Function Curl wrapper
     *
     * This class template represents the curl of another vector field.
     *
     * \tparam Function_
     * The function whose curl is to be wrapped.
     *
     * \author Peter Zajac
     */
    template<typename Function_>
    class Curl :
      public Analytic::Function
    {
    public:
      // ensure that the input function is scalar
      static_assert(Function_::ImageType::is_vector, "only vector fields are supported; use ScalarCurl for scalar functions");

      /// our domain dimension is the same as the input function's
      static constexpr int domain_dim = Function_::domain_dim;

      /// this is a vector-valued function
      typedef Image::Vector<domain_dim> ImageType;

      static constexpr bool can_value = Function_::can_grad;
      static constexpr bool can_grad = Function_::can_hess;
      static constexpr bool can_hess = false;

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;
        typedef typename Traits_::HessianType HessianType;

      private:
        typedef EvalTraits<DataType, Function_> FuncEvalTraits;
        /// our original function evaluator
        typename Function_::template Evaluator<FuncEvalTraits> _func_eval;

        /// 2D vector curl operator
        template<typename T_, int s_, int sm_, int sn_>
        static void compute(Tiny::Vector<T_, 2, s_>& curl, const Tiny::Matrix<T_, 2, 2, sm_, sn_>& grad)
        {
          curl[0] = -grad[1][0];
          curl[1] = +grad[0][1];
        }

        template<typename T_, int n_, int sm1_, int sn1_, int sl_, int sm_, int sn_>
        static void compute(Tiny::Matrix<T_, 2, n_, sm1_, sn1_>& curl, const Tiny::Tensor3<T_, 2, 2, n_, sl_, sm_, sn_>& grad)
        {
          curl[0] = -grad[1][0];
          curl[1] = +grad[0][1];
        }

        /// 3D vector curl operator
        template<typename T_, int s_, int sm_, int sn_>
        static void compute(Tiny::Vector<T_, 3, s_>& curl, const Tiny::Matrix<T_, 3, 3, sm_, sn_>& grad)
        {
          curl[0] = grad[1][2] - grad[2][1];
          curl[1] = grad[2][0] - grad[0][2];
          curl[2] = grad[0][1] - grad[1][0];
        }

        template<typename T_, int n_, int sm1_, int sn1_, int sl_, int sm_, int sn_>
        static void compute(Tiny::Matrix<T_, 3, n_, sm1_, sn1_>& curl, const Tiny::Tensor3<T_, 3, 3, n_, sl_, sm_, sn_>& grad)
        {
          curl[0] = grad[1][2] - grad[2][1];
          curl[1] = grad[2][0] - grad[0][2];
          curl[2] = grad[0][1] - grad[1][0];
        }

      public:
        explicit Evaluator(const Curl& function) :
          _func_eval(function._function)
        {
        }

        void value(ValueType& val, const PointType& point)
        {
          typename FuncEvalTraits::GradientType grad;
          _func_eval.gradient(grad, point);
          compute(val, grad);
        }

        void gradient(GradientType& grad, const PointType& point)
        {
          typename FuncEvalTraits::HessianType hess;
          _func_eval.hessian(hess, point);
          compute(grad, hess);
        }
      };

    private:
      const Function_& _function;

    public:
      explicit Curl(const Function_& function) :
        _function(function)
      {
      }
    };

    /**
     * \brief Analytic Scalar Function Curl wrapper
     *
     * This class template represents the curl of a scalar function.
     *
     * \tparam Function_
     * The function whose curl is to be wrapped.
     *
     * \author Peter Zajac
     */
    template<typename Function_>
    class ScalarCurl :
      public Analytic::Function
    {
    public:
      // ensure that the input function is scalar
      static_assert(Function_::ImageType::is_scalar, "only scalar functions are supported; use Curl for vector fields");

      /// our domain dimension is the same as the input function's
      static constexpr int domain_dim = Function_::domain_dim;

      /// this is a vector-valued function
      typedef Image::Vector<domain_dim> ImageType;

      static constexpr bool can_value = Function_::can_grad;
      static constexpr bool can_grad = Function_::can_hess;
      static constexpr bool can_hess = false;

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;
        typedef typename Traits_::HessianType HessianType;

      private:
        typedef EvalTraits<DataType, Function_> FuncEvalTraits;
        /// our original function evaluator
        typename Function_::template Evaluator<FuncEvalTraits> _func_eval;

        /// 2D vector curl operator
        template<typename T_, int s_, int sn_>
        static void compute(Tiny::Vector<T_, 2, s_>& curl, const Tiny::Vector<T_, 2, sn_>& grad)
        {
          curl[0] = -grad[1];
          curl[1] = +grad[0];
        }

        template<typename T_, int sa_, int sb_, int sm_, int sn_>
        static void compute(Tiny::Matrix<T_, 2, 2, sa_, sb_>& curl, const Tiny::Matrix<T_, 2, 2, sm_, sn_>& grad)
        {
          curl[0] = -grad[1];
          curl[1] = +grad[0];
        }

      public:
        explicit Evaluator(const ScalarCurl& function) :
          _func_eval(function._function)
        {
        }

        void value(ValueType& val, const PointType& point)
        {
          typename FuncEvalTraits::GradientType grad;
          _func_eval.gradient(grad, point);
          compute(val, grad);
        }

        void gradient(GradientType& grad, const PointType& point)
        {
          typename FuncEvalTraits::HessianType hess;
          _func_eval.hessian(hess, point);
          compute(grad, hess);
        }
      };

    private:
      const Function_& _function;

    public:
      explicit ScalarCurl(const Function_& function) :
        _function(function)
      {
      }
    };
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_WRAPPERS_HPP
