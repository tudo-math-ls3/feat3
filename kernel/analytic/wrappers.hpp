// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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

        ValueType value(const PointType& point)
        {
          return _func_eval.gradient(point);
        }

        GradientType gradient(const PointType& point)
        {
          return _func_eval.hessian(point);
        }
      };

    private:
      const Function_& _function;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the function whose gradient is to be wrapped.
       */
      explicit Gradient(const Function_& function) :
        _function(function)
      {
      }
    }; // class Gradient<...>

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

        ValueType value(const PointType& point)
        {
          return _func_eval.gradient(point).trace();
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
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the function whose divergence is to be wrapped.
       */
      explicit Divergence(const Function_& function) :
        _function(function)
      {
      }
    }; // class Divergence<...>

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
        static Tiny::Vector<T_, 2, s_> compute(const Tiny::Matrix<T_, 2, 2, sm_, sn_>& grad)
        {
          Tiny::Vector<T_, 2, s_> curl;
          curl[0] = -grad[1][0];
          curl[1] = +grad[0][1];
          return curl;
        }

        template<typename T_, int n_, int sm1_, int sn1_, int sl_, int sm_, int sn_>
        static Tiny::Matrix<T_, 2, n_, sm1_, sn1_> compute(const Tiny::Tensor3<T_, 2, 2, n_, sl_, sm_, sn_>& grad)
        {
          Tiny::Matrix<T_, 2, n_, sm1_, sn1_> curl;
          curl[0] = -grad[1][0];
          curl[1] = +grad[0][1];
          return curl;
        }

        /// 3D vector curl operator
        template<typename T_, int s_, int sm_, int sn_>
        static Tiny::Vector<T_, 3, s_> compute(const Tiny::Matrix<T_, 3, 3, sm_, sn_>& grad)
        {
          Tiny::Vector<T_, 3, s_> curl;
          curl[0] = grad[1][2] - grad[2][1];
          curl[1] = grad[2][0] - grad[0][2];
          curl[2] = grad[0][1] - grad[1][0];
          return curl;
        }

        template<typename T_, int n_, int sm1_, int sn1_, int sl_, int sm_, int sn_>
        static Tiny::Matrix<T_, 3, n_, sm1_, sn1_> compute(const Tiny::Tensor3<T_, 3, 3, n_, sl_, sm_, sn_>& grad)
        {
          Tiny::Matrix<T_, 3, n_, sm1_, sn1_> curl;
          curl[0] = grad[1][2] - grad[2][1];
          curl[1] = grad[2][0] - grad[0][2];
          curl[2] = grad[0][1] - grad[1][0];
          return curl;
        }

      public:
        explicit Evaluator(const Curl& function) :
          _func_eval(function._function)
        {
        }

        ValueType value(const PointType& point)
        {
          return compute(_func_eval.gradient(point));
        }

        GradientType gradient(const PointType& point)
        {
          return compute(_func_eval.hessian(point));
        }
      };

    private:
      const Function_& _function;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the function whose curl is to be wrapped.
       */
      explicit Curl(const Function_& function) :
        _function(function)
      {
      }
    }; // class Curl<...>

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

        /// 2D scalar curl operator
        template<typename T_, int s_, int sn_>
        static Tiny::Vector<T_, 2, s_> compute(const Tiny::Vector<T_, 2, sn_>& grad)
        {
          Tiny::Vector<T_, 2, s_> curl;
          curl[0] = -grad[1];
          curl[1] = +grad[0];
          return curl;
        }

        template<typename T_, int sa_, int sb_, int sm_, int sn_>
        static Tiny::Matrix<T_, 2, 2, sa_, sb_> compute(const Tiny::Matrix<T_, 2, 2, sm_, sn_>& grad)
        {
          Tiny::Matrix<T_, 2, 2, sa_, sb_> curl;
          curl[0] = -grad[1];
          curl[1] = +grad[0];
          return curl;
        }

        /// 3D scalar curl operator
        template<typename T_, int s_, int sn_>
        static Tiny::Vector<T_, 3, s_> compute(const Tiny::Vector<T_, 3, sn_>& grad)
        {
          Tiny::Vector<T_, 3, s_> curl;
          curl[0] = grad[1] - grad[2];
          curl[1] = grad[2] - grad[0];
          curl[2] = grad[0] - grad[1];
          return curl;
        }

        template<typename T_, int sa_, int sb_, int sm_, int sn_>
        static Tiny::Matrix<T_, 3, 3, sa_, sb_> compute(const Tiny::Matrix<T_, 3, 3, sm_, sn_>& grad)
        {
          Tiny::Matrix<T_, 3, 3, sa_, sb_> curl;
          curl[0] = grad[1] - grad[2];
          curl[1] = grad[2] - grad[0];
          curl[2] = grad[0] - grad[1];
          return curl;
        }

      public:
        explicit Evaluator(const ScalarCurl& function) :
          _func_eval(function._function)
        {
        }

        ValueType value(const PointType& point)
        {
          return compute(_func_eval.gradient(point));
        }

        GradientType gradient(const PointType& point)
        {
          return compute(_func_eval.hessian(point));
        }
      };

    private:
      const Function_& _function;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the function whose scalar is to be wrapped.
       */
      explicit ScalarCurl(const Function_& function) :
        _function(function)
      {
      }
    }; // class ScalarCurl<...>

    /**
     *  \brief This class is a wrapper transforming a polar-basis function to a euclidean base one
     *
     *  \tparam Function_ The function object to wrap
     *  \tparam pos_range Should we use angles in \f$ [0, 2\pi] \f$ ?
     *
     *  This wrapper handles the transformation of a 2D function defined in (shifted and rotated) polar coordinates,
     *  i.e. function defined simply in \f$ (r,\varphi) \f$
     *  to the standard euclidean basis \f$ x,y \f$, where the polar coordinates can be expressed by a orthogonally transformed euclidean basis
     *  \f$(e_x', e_y')\f$ with \f$ x' = r \cos(\varphi), y' = r \sin(\varphi) \f$. To be more precisly:
     *  This Wrapper expresses values, gradients and hessians (if available) in the standard euclidean basis e_x, e_y
     *  while the actual function can be defined in an arbitrarly rotated and shifted polar-coordinate space.
     *  For the affin-linear transformation fa(x) = f(Ax+c), we simply use the chain rule and arrive at
     *          grad fa(x) = A^T grad f(z),
     *          hess fa(x) = A^T hess f(z) A,
     *   where z = Ax+c.
     *  The influence of the change to curveilinear coordinates, that is hidden inside of f is a bit more complicated
     *  but can be arrived at using the simple connection of
     *  \f$ du/dx = \frac{\partial u}{\partial r} \frac{\partial r}{\partial x} + \frac{\partial u}{\partial \varphi} \frac{\partial \varphi}{\partial x} \f$ .
     *
     *  \author Maximilian Esser
     */
    template<typename Function_, bool pos_range_ = true>
    class PolarCoordinate :
      public Analytic::Function
    {
    public:
      // ensure that the input function is scalar
      static_assert(Function_::ImageType::is_scalar, "For now, we only allow scalar functions... although this should be doable without any major problems.");

      /// our domain dimension is the same as the input function's
      static constexpr int domain_dim = Function_::domain_dim;
      static constexpr bool pos_range = pos_range_;

      static_assert(domain_dim == 2, "PolarCoordinates are only available for 2 dimenions!");

      /// this is a vector-valued function
      typedef typename Function_::ImageType ImageType;
      typedef typename Function_::DataType DataType;

      static constexpr bool can_value = Function_::can_value;
      static constexpr bool can_grad = Function_::can_grad;
      static constexpr bool can_hess = Function_::can_hess;

      //our rotation matrix defined as x = rot * x' <- ie. the inverse rotation to the theta given in the ctors
      Tiny::Matrix<DataType, 2, 2> _rotation;
      Tiny::Vector<DataType, 2> _shift;

      PolarCoordinate(const Function_& func) :
      _rotation(),
      _shift(Tiny::Vector<DataType, 2>::null()),
      _function(func)
      {
        _rotation.set_identity();
      }

      /// theta rot angle, so that x' = rot(theta) * (x - shift)
      /// in case you actually want to rotate the coordinate system, keep in mind,
      /// that this is done by using -theta...
      PolarCoordinate(const Function_& func, DataType theta) :
      _rotation(),
      _shift(Tiny::Vector<DataType, 2>::null()),
      _function(func)
      {
        _rotation.set_rotation_2d(-theta);
      }

      /// theta rot angle, so that x' = rot(theta) * x
      PolarCoordinate(const Function_& func, Tiny::Vector<DataType, 2> shift, DataType theta = DataType(0)) :
      _rotation(),
      _shift(std::move(shift)),
      _function(func)
      {
        _rotation.set_rotation_2d(-theta);
      }

      /// provide shift vector and the first x' basis vector
      PolarCoordinate(const Function_& func, Tiny::Vector<DataType, 2> shift, Tiny::Vector<DataType, 2> x_base) :
      _rotation(),
      _shift(std::move(shift)),
      _function(func)
      {
        // theta to rotate (1,0) into x_base
        DataType rot_theta = Math::calc_opening_angle(DataType(1), DataType(0), x_base[0], x_base[1]);
        _rotation.set_rotation_2d(-rot_theta);

      }

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
        const Tiny::Matrix<DataType, 2, 2>& _rotation;
        Tiny::Matrix<DataType, 2, 2> _rotation_t;
        const Tiny::Vector<DataType, 2>& _shift;


      PointType transform(const PointType& point)
      {
        PointType pcoords;
        //shift and rotate
        auto p_n = _rotation * (point - _shift);
        //and now calculate the actual polar coordinates
        pcoords[0] = p_n.norm_euclid();
        pcoords[1] = (pcoords[0] < DataType(5)*Math::eps<DataType>()) ? DataType(0) : Math::acos(p_n[0]/pcoords[0]);
        //annoying fix in case of float precision...
        if constexpr(pos_range_)
          pcoords[1] = ((p_n[1] >= DataType(0)) || (pcoords[1] < DataType(3)*Math::eps<DataType>())) ? pcoords[1] : DataType(2)*Math::pi<DataType>()-pcoords[1];
        else
          pcoords[1] = (p_n[1] < DataType(0)) ? -pcoords[1] : pcoords[1];
        return pcoords;
      }

      GradientType compute(const GradientType& grad_p, const PointType& pcoords)
      {
        GradientType grad;
        grad[0] = grad_p[0] * Math::cos(pcoords[1]);
        grad[1] = grad_p[0] * Math::sin(pcoords[1]);
        if(pcoords[0] > DataType(5) * Math::eps<DataType>())
        {
          grad[0] -= grad_p[1] * Math::sin(pcoords[1]) / pcoords[0];
          grad[1] += grad_p[1] * Math::cos(pcoords[1]) / pcoords[0];
        }
        return _rotation_t * grad;
      }

      HessianType compute(const HessianType& hess_p, const GradientType& grad_p, const PointType& pcoords)
      {
        HessianType hess;
        const DataType valsin = Math::sin(pcoords[1]);
        const DataType valcos = Math::cos(pcoords[1]);
        hess[0][0] = hess_p[0][0] * Math::sqr(valcos);
        hess[1][1] = hess_p[0][0] * Math::sqr(valsin);
        hess[0][1] = hess_p[0][0] * valsin * valcos;
        if(pcoords[0] > DataType(5) * Math::sqrt(Math::eps<DataType>()))
        {
          const DataType r_re = DataType(1) / pcoords[0];
          hess[0][0] += grad_p[0] * Math::sqr(valsin) * r_re + DataType(2) * grad_p[1] * valsin * valcos * Math::sqr(r_re)
                          - DataType(2) * hess_p[0][1] * valsin * valcos * r_re + hess_p[1][1] * Math::sqr(valsin * r_re);
          hess[1][1] += grad_p[0] * Math::sqr(valcos) * r_re - DataType(2) * grad_p[1] * valsin * valcos * Math::sqr(r_re)
                          + DataType(2) * hess_p[0][1] * valsin * valcos * r_re + hess_p[1][1] * Math::sqr(valcos * r_re);
          hess[0][1] += - grad_p[0] * valsin * valcos * r_re + grad_p[1] * (Math::sqr(valsin) - Math::sqr(valcos)) * Math::sqr(r_re)
                          + hess_p[0][1] * (Math::sqr(valcos) - Math::sqr(valsin)) * r_re - hess_p[1][1] * valsin * valcos * Math::sqr(r_re);

        }
        hess[1][0] = hess[0][1];
        return _rotation_t * hess * _rotation;

      }

      public:
        explicit Evaluator(const PolarCoordinate& function) :
          _func_eval(function._function),
          _rotation(function._rotation),
          _shift(function._shift)
        {
          _rotation_t.set_transpose(_rotation);
        }

        ValueType value(const PointType& point)
        {
          return _func_eval.value(this->transform(point));
        }

        GradientType gradient(const PointType& point)
        {
          auto p = this->transform(point);
          return compute(_func_eval.gradient(p), p);
        }

        HessianType hessian(const PointType& point)
        {
          auto p = this->transform(point);
          return compute(_func_eval.hessian(p), _func_eval.gradient(p), p);
        }
      };

      private:
      const Function_& _function;
    }; // class PolarCoordinate<...>

  } // namespace Analytic
} // namespace FEAT
