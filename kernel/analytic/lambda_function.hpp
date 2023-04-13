// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_LAMBDA_FUNCTION_HPP
#define KERNEL_ANALYTIC_LAMBDA_FUNCTION_HPP 1

#include <kernel/analytic/function.hpp>
#include <kernel/util/exception.hpp>

namespace FEAT
{
  namespace Analytic
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Lambda Expression container helper class
       *
       * This helper class is used by the LambdaSet class templates to store the actual lambda
       * expressions supplied by the user -- if the lambda is not void, in which case a specialization
       * of this class is used.
       *
       * \tparam F_
       * The lamdba expression type as supplied by the user or \c void, if no lambda expression was given
       *
       * \author Peter Zajac
       */
      template<typename F_>
      struct LambdaHelper
      {
        /// the lambda expression type; not void
        typedef F_ Type;
        /// we have an expression that's not void (otherwise the specialization would have been chosen)
        static constexpr bool have = true;
        /// the actual lambda expression functor
        F_ f;
        /// a handy constructor
        explicit LambdaHelper(F_&& f_) : f(std::forward<F_>(f_)) {}
      };

      /**
       * \brief Specialization for when a lambda expression is not given aka. void
       */
      template<>
      struct LambdaHelper<void>
      {
        /// no lambda to be found here
        typedef void Type;
        /// we're out of lambda
        static constexpr bool have = false;
        /// standard constructor
        LambdaHelper() {}
      };

      /**
       * \brief Helper class for the evaluation of a gradient by using explicitly given lambda expressions
       *
       * This helper class is used to evaluate a gradient by simply evaluating the lambda expressions
       * that the user has supplied for the partial derivatives.
       *
       * \tparam LambdaSet_
       * The lamdba set that is stored in the Lambda function's evaluator
       *
       * \tparam dim_
       * The domain/image dimension
       *
       * \tparam DT_
       * The fundamental floating point data type
       *
       * \tpatam VT_
       * The function value type:
       * - equal to DT_ for scalar functions
       * - equal to Tiny::Vector<DT_,dim_> for vector fields
       *
       * \tparam GT_
       * The function gradient type:
       * - equal to Tiny::Vector<DT_,dim_> for scalar functions
       * - equal to Tiny::Matrix<DT_,dim_,dim_> for vector fields
       *
       * \tparam have_grad_
       * Set to \c true, if the user has suppled all first-order partial derivatives required for the dimension
       * as lambda expressions, otherwise \c false
       *
       * \author Peter Zajac
       */
      template<typename LambdaSet_, int dim_, typename DT_, typename VT_, typename GT_, bool have_grad_>
      class LambdaGradHelper
      {
      public:
        /// the lambda set of the evaluator
        LambdaSet_& lambda_set;

        explicit LambdaGradHelper(LambdaSet_& lambda_set_, DT_, int) :
          lambda_set(lambda_set_)
        {
        }

        /// overload for scalar 1D gradients
        void eval_grad(Tiny::Vector<DT_, 1>& grad, const Tiny::Vector<DT_, 1>& p)
        {
          grad[0] = lambda_set.dx.f(p[0]);
        }

        /// overload for scalar 2D gradients
        void eval_grad(Tiny::Vector<DT_, 2>& grad, const Tiny::Vector<DT_, 2>& p)
        {
          grad[0] = lambda_set.dx.f(p[0], p[1]);
          grad[1] = lambda_set.dy.f(p[0], p[1]);
        }

        /// overload for scalar 3D gradients
        void eval_grad(Tiny::Vector<DT_, 3>& grad, const Tiny::Vector<DT_, 3>& p)
        {
          grad[0] = lambda_set.dx.f(p[0], p[1], p[2]);
          grad[1] = lambda_set.dy.f(p[0], p[1], p[2]);
          grad[2] = lambda_set.dz.f(p[0], p[1], p[2]);
        }

        /// overload for vector 2D gradients
        void eval_grad(Tiny::Matrix<DT_, 2, 2>& grad, const Tiny::Vector<DT_, 2>& p)
        {
          grad[0][0] = lambda_set.dx1.f(p[0], p[1]);
          grad[0][1] = lambda_set.dy1.f(p[0], p[1]);
          grad[1][0] = lambda_set.dx2.f(p[0], p[1]);
          grad[1][1] = lambda_set.dy2.f(p[0], p[1]);
        }

        /// overload for vector 3D gradients
        void eval_grad(Tiny::Matrix<DT_, 3, 3>& grad, const Tiny::Vector<DT_, 3>& p)
        {
          grad[0][0] = lambda_set.dx1.f(p[0], p[1], p[2]);
          grad[0][1] = lambda_set.dy1.f(p[0], p[1], p[2]);
          grad[0][2] = lambda_set.dz1.f(p[0], p[1], p[2]);
          grad[1][0] = lambda_set.dx2.f(p[0], p[1], p[2]);
          grad[1][1] = lambda_set.dy2.f(p[0], p[1], p[2]);
          grad[1][2] = lambda_set.dz2.f(p[0], p[1], p[2]);
          grad[2][0] = lambda_set.dx3.f(p[0], p[1], p[2]);
          grad[2][1] = lambda_set.dy3.f(p[0], p[1], p[2]);
          grad[2][2] = lambda_set.dz3.f(p[0], p[1], p[2]);
        }
      }; // class LambdaGradHelper<...> (explicit formulae)

      /**
       * \brief Specialization of the gradient evaluation helper class for Richardson extrapolation
       *
       * This specialization is chosen if the user did not supply all (or any) first-order partial derivatives as
       * lambda expressions and therefore the evaluation of the gradient has to be performed by Richardson extrapolation
       * of second-order central difference quotients, which is what this helper class is actually responsible for.
       *
       * \author Peter Zajac
       */
      template<typename LambdaSet_, int dim_, typename DT_, typename VT_, typename GT_>
      class LambdaGradHelper<LambdaSet_, dim_, DT_, VT_, GT_, false>
      {
      public:
        /// the lambda set of the evaluator
        LambdaSet_& lambda_set;
        /// the initial h for the Richardson extrapolation; usually 0.01
        const DT_ initial_h;
        /// the maximum number of Richardson extrapolation steps; usually 10
        const std::size_t max_steps;

        explicit LambdaGradHelper(LambdaSet_& lambda_set_, DT_ initial_h_, int max_steps_) :
          lambda_set(lambda_set_),
          initial_h(initial_h_),
          max_steps(std::size_t(max_steps_))
        {
        }

        /// extrapolation method for gradients (all dimensions)
        void eval_grad(GT_& grad, const Tiny::Vector<DT_, dim_>& point)
        {
          // first, create a mutable copy of our point
          Tiny::Vector<DT_, dim_> v(point);

          /// our gradient extrapolation table
          std::vector<GT_> gradex(max_steps);

          // next, choose the initial h
          DT_ h(initial_h);

          // evaluate gradient
          eval_grad_diff_quot(gradex[0], v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(gradex.size()-1);
          DT_ def(Math::huge<DT_>());
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next gradient
            eval_grad_diff_quot(gradex[i+1], v, h *= DT_(0.5));

            // initialize scaling factor
            DT_ q = DT_(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DT_(4);
              (gradex[k] -= q*gradex[k+1]) *= DT_(1) / (DT_(1) - q);
            }

            // compute and check our defect
            DT_ d = def_norm_sqr(gradex[1], gradex[0]);
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              break;
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DT_(4);
            (gradex[0] -= q*gradex[1]) *= DT_(1) / (DT_(1) - q);
          }

          // return our extrapolated gradient
          grad = gradex.front();
        }

        /// squared defect norm for scalar gradients
        template<int n_, int s_>
        static DT_ def_norm_sqr(const Tiny::Vector<DT_, n_, s_>& x, const Tiny::Vector<DT_, n_, s_>& y)
        {
          return (x - y).norm_euclid_sqr();
        }

        /// squared defect norm for vector gradients
        template<int m_, int n_, int sm_, int sn_>
        static DT_ def_norm_sqr(
          const Tiny::Matrix<DT_, m_, n_, sm_, sn_>& x,
          const Tiny::Matrix<DT_, m_, n_, sm_, sn_>& y)
        {
          return (x - y).norm_hessian_sqr();
        }

        // compute scalar difference quotient
        static void set_diff_quot(Tiny::Vector<DT_, dim_>& x, int i, const VT_& vr, const VT_& vl, const DT_ denom)
        {
          x[i] = denom * (vr - vl);
        }

        // compute vector difference quotient
        static void set_diff_quot(Tiny::Matrix<DT_, dim_, dim_>& x, int i, const VT_& vr, const VT_& vl, const DT_ denom)
        {
          for(int k(0); k < dim_; ++k)
            x[k][i] = denom * (vr[k] - vl[k]);
        }

        /// evaluates a central difference quotient for the first derivative (all dimensions)
        void eval_grad_diff_quot(GT_& x, Tiny::Vector<DT_, dim_>& v, const DT_ h)
        {
          // difference quotient denominator
          const DT_ denom = DT_(1) / (DT_(2) * h);
          VT_ vr, vl;

          // loop over all dimensions
          for(int i(0); i < dim_; ++i)
          {
            // backup coordinate i
            const DT_ vi(v[i]);

            // evaluate f(v + h*e_i)
            v[i] = vi + h;
            vr = lambda_set.eval_value(v);

            // evaluate f(v - h*e_i)
            v[i] = vi - h;
            vl = lambda_set.eval_value(v);

            // compute difference quotient
            set_diff_quot(x, i, vr, vl, denom);

            // restore coord
            v[i] = vi;
          }
        }
      }; // class LambdaGradHelper<...> (Richardson extrapolation)

      /**
       * \brief Helper class for the evaluation of Hessians by using explicitly given lambda expressions
       *
       * This helper class is used to evaluate Hessians by simply evaluating the lambda expressions
       * that the user has supplied for the partial derivatives.
       *
       * \tparam LambdaSet_
       * The lamdba set that is stored in the Lambda function's evaluator
       *
       * \tparam dim_
       * The domain/image dimension
       *
       * \tparam DT_
       * The fundamental floating point data type
       *
       * \tpatam VT_
       * The function value type:
       * - equal to DT_ for scalar functions
       * - equal to Tiny::Vector<DT_,dim_> for vector fields
       *
       * \tparam HT_
       * The function Hessian type:
       * - equal to Tiny::Matrix<DT_,dim_,dim_> for scalar functions
       * - equal to Tiny::Tensor3<DT_,dim_,dim_,dim_> for vector fields
       *
       * \tparam have_hess_
       * Set to \c true, if the user has suppled all second-order partial derivatives required for the dimension
       * as lambda expressions, otherwise \c false
       *
       * \author Peter Zajac
       */
      template<typename LambdaSet_, int dim_, typename DT_, typename VT_, typename HT_, bool have_hess_>
      class LambdaHessHelper
      {
      public:
        /// the lambda set of the evaluator
        LambdaSet_& lambda_set;

        explicit LambdaHessHelper(LambdaSet_& lambda_set_, DT_, int) :
          lambda_set(lambda_set_)
        {
        }

        /// overload for scalar 1D hessians
        void eval_hess(Tiny::Matrix<DT_, 1, 1>& hess, const Tiny::Vector<DT_, 1>& p)
        {
          hess[0][0] = lambda_set.dxx.f(p[0]);
        }

        /// overload for scalar 2D hessians
        void eval_hess(Tiny::Matrix<DT_, 2, 2>& hess, const Tiny::Vector<DT_, 2>& p)
        {
          hess[0][0] = lambda_set.dxx.f(p[0], p[1]);
          hess[1][1] = lambda_set.dyy.f(p[0], p[1]);
          hess[0][1] = hess[1][0] = lambda_set.dxy.f(p[0], p[1]);
        }

        /// overload for scalar 3D hessians
        void eval_hess(Tiny::Matrix<DT_, 3, 3>& hess, const Tiny::Vector<DT_, 3>& p)
        {
          hess[0][0] = lambda_set.dxx.f(p[0], p[1], p[2]);
          hess[1][1] = lambda_set.dyy.f(p[0], p[1], p[2]);
          hess[2][2] = lambda_set.dzz.f(p[0], p[1], p[2]);
          hess[0][1] = hess[1][0] = lambda_set.dxy.f(p[0], p[1], p[2]);
          hess[1][2] = hess[2][1] = lambda_set.dyz.f(p[0], p[1], p[2]);
          hess[2][0] = hess[0][2] = lambda_set.dzx.f(p[0], p[1], p[2]);
        }

        /// overload for vector 2D hessians
        void eval_hess(Tiny::Tensor3<DT_, 2, 2, 2>& hess, const Tiny::Vector<DT_, 2>& p)
        {
          hess[0][0][0] = lambda_set.dxx1.f(p[0], p[1]);
          hess[0][1][1] = lambda_set.dyy1.f(p[0], p[1]);
          hess[0][0][1] = hess[0][1][0] = lambda_set.dxy1.f(p[0], p[1]);
          hess[1][0][0] = lambda_set.dxx2.f(p[0], p[1]);
          hess[1][1][1] = lambda_set.dyy2.f(p[0], p[1]);
          hess[1][0][1] = hess[1][1][0] = lambda_set.dxy2.f(p[0], p[1]);
        }

        /// overload for vector 3D hessians
        void eval_hess(Tiny::Tensor3<DT_, 3, 3, 3>& hess, const Tiny::Vector<DT_, 3>& p)
        {
          hess[0][0][0] = lambda_set.dxx1.f(p[0], p[1], p[2]);
          hess[0][1][1] = lambda_set.dyy1.f(p[0], p[1], p[2]);
          hess[0][2][2] = lambda_set.dzz1.f(p[0], p[1], p[2]);
          hess[0][0][1] = hess[0][1][0] = lambda_set.dxy1.f(p[0], p[1], p[2]);
          hess[0][1][2] = hess[0][2][1] = lambda_set.dyz1.f(p[0], p[1], p[2]);
          hess[0][2][0] = hess[0][0][2] = lambda_set.dzx1.f(p[0], p[1], p[2]);
          hess[1][0][0] = lambda_set.dxx2.f(p[0], p[1], p[2]);
          hess[1][1][1] = lambda_set.dyy2.f(p[0], p[1], p[2]);
          hess[1][2][2] = lambda_set.dzz2.f(p[0], p[1], p[2]);
          hess[1][0][1] = hess[1][1][0] = lambda_set.dxy2.f(p[0], p[1], p[2]);
          hess[1][1][2] = hess[1][2][1] = lambda_set.dyz2.f(p[0], p[1], p[2]);
          hess[1][2][0] = hess[1][0][2] = lambda_set.dzx2.f(p[0], p[1], p[2]);
          hess[2][0][0] = lambda_set.dxx3.f(p[0], p[1], p[2]);
          hess[2][1][1] = lambda_set.dyy3.f(p[0], p[1], p[2]);
          hess[2][2][2] = lambda_set.dzz3.f(p[0], p[1], p[2]);
          hess[2][0][1] = hess[2][1][0] = lambda_set.dxy3.f(p[0], p[1], p[2]);
          hess[2][1][2] = hess[2][2][1] = lambda_set.dyz3.f(p[0], p[1], p[2]);
          hess[2][2][0] = hess[2][0][2] = lambda_set.dzx3.f(p[0], p[1], p[2]);
        }
      }; // class LambdaHessHelper<...> (explicit 2D formula)

      /**
       * \brief Specialization of the Hessian evaluation helper class for Richardson extrapolation
       *
       * This specialization is chosen if the user did not supply all (or any) second-order partial derivatives as
       * lambda expressions and therefore the evaluation of the Hessian has to be performed by Richardson extrapolation
       * of second-order central difference quotients, which is what this helper class is actually responsible for.
       *
       * \author Peter Zajac
       */
      template<typename LambdaSet_, int dim_, typename DT_, typename VT_, typename HT_>
      class LambdaHessHelper<LambdaSet_, dim_, DT_, VT_, HT_, false>
      {
      public:
        LambdaSet_& lambda_set;
        const DT_ initial_h;
        const std::size_t max_steps;

        explicit LambdaHessHelper(LambdaSet_& lambda_set_, DT_ initial_h_, int max_steps_) :
          lambda_set(lambda_set_),
          initial_h(initial_h_),
          max_steps(std::size_t(max_steps_))
        {
        }

        void eval_hess(HT_& hess, const Tiny::Vector<DT_, dim_>& point)
        {
          // first, create a mutable copy of our point
          Tiny::Vector<DT_, dim_> v(point);

          /// our gradient extrapolation table
          std::vector<HT_> hessex(max_steps);

          // next, choose the initial h
          DT_ h(initial_h);

          // evaluate hessian
          eval_hess_diff_quot(hessex[0], v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(hessex.size()-1);
          DT_ def(Math::huge<DT_>());
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next hessian
            eval_hess_diff_quot(hessex[i+1], v, h *= DT_(0.5));

            // initialize scaling factor
            DT_ q = DT_(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DT_(4);
              (hessex[k] -= q*hessex[k+1]) *= DT_(1) / (DT_(1) - q);
            }

            // compute and check our defect
            DT_ d = def_norm_sqr(hessex[1], hessex[0]);
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              break;
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DT_(4);
            (hessex[0] -= q*hessex[1]) *= DT_(1) / (DT_(1) - q);
          }

          // return our extrapolated hessian
          hess = hessex.front();
        }

        /// squared defect norm for scalar hessians
        template<int m_, int n_, int sm_, int sn_>
        static DT_ def_norm_sqr(
          const Tiny::Matrix<DT_, m_, n_, sm_, sn_>& x,
          const Tiny::Matrix<DT_, m_, n_, sm_, sn_>& y)
        {
          return (x - y).norm_hessian_sqr();
        }

        /// squared defect norm for vector hessians
        template<int l_, int m_, int n_, int sl_, int sm_, int sn_>
        static DT_ def_norm_sqr(
          const Tiny::Tensor3<DT_, l_, m_, n_, sl_, sm_, sn_>& x,
          const Tiny::Tensor3<DT_, l_, m_, n_, sl_, sm_, sn_>& y)
        {
          DT_ r(DT_(0));
          for(int i(0); i < l_; ++i)
          {
            for(int j(0); j < m_; ++j)
            {
              for(int k(0); k < n_; ++k)
              {
                r += Math::sqr(x(i,j,k) - y(i,j,k));
              }
            }
          }
          return r;
        }

        // compute scalar difference quotient for i == j
        static void set_diff_quot(Tiny::Matrix<DT_, dim_, dim_>& x, int i, const VT_& vr, const VT_& vl,  const VT_& vc, const DT_ denom)
        {
          x[i][i] = denom * (vr + vl - DT_(2)*vc);
        }

        // compute scalar difference quotient for i != j
        static void set_diff_quot(Tiny::Matrix<DT_, dim_, dim_>& x, int i, int j, const VT_& vne, const VT_& vsw, const VT_& vnw, const VT_& vse, const DT_ denom)
        {
          x[i][j] = x[j][i] = denom * ((vne + vsw) - (vnw + vse));
        }

        // compute vector difference quotient for i == j
        static void set_diff_quot(Tiny::Tensor3<DT_, dim_, dim_, dim_>& x, int i, const VT_& vr, const VT_& vl,  const VT_& vc, const DT_ denom)
        {
          for(int k(0); k < dim_; ++k)
            x[k][i][i] = denom * (vr[k] + vl[k] - DT_(2)*vc[k]);
        }

        // compute vector difference quotient for i != j
        static void set_diff_quot(Tiny::Tensor3<DT_, dim_, dim_, dim_>& x, int i, int j, const VT_& vne, const VT_& vsw, const VT_& vnw, const VT_& vse, const DT_ denom)
        {
          for(int k(0); k < dim_; ++k)
            x[k][i][j] = x[k][j][i] = denom * ((vne[k] + vsw[k]) - (vnw[k] + vse[k]));
        }

        /// evaluates the second-order difference quotients
        void eval_hess_diff_quot(HT_& x, Tiny::Vector<DT_, dim_>& v, const DT_ h)
        {
          // difference quotient denominators
          const DT_ denom1 = DT_(1) / (h * h);
          const DT_ denom2 = DT_(1) / (DT_(4) * h * h);
          VT_ vc, vr, vl, vne, vnw, vse, vsw;

          // evaluate at point
          vc = lambda_set.eval_value(v);

          // loop over all dimensions
          for(int i(0); i < dim_; ++i)
          {
            // backup coord
            const DT_ vi(v[i]);

            // eval f(x + h*e_i)
            v[i] = vi + h;
            vr = lambda_set.eval_value(v);

            // eval f(x-h)
            v[i] = vi - h;
            vl = lambda_set.eval_value(v);

            // compute difference quotient
            set_diff_quot(x, i, vr, vl, vc, denom1);

            // now the mixed derivatives
            for(int j(i+1); j < dim_; ++j)
            {
              // backup coord
              const DT_ vj(v[j]);

              // we need four points here:
              // north-east: f(v + h*e_i + h*e_j)
              v[i] = vi + h;
              v[j] = vj + h;
              vne = lambda_set.eval_value(v);

              // north-west: f(v - h*e_i + h*e_j)
              v[i] = vi - h;
              v[j] = vj + h;
              vnw = lambda_set.eval_value(v);

              // south-east: f(v + h*e_i - h*e_j)
              v[i] = vi + h;
              v[j] = vj - h;
              vse = lambda_set.eval_value(v);

              // south-west: f(v - h*e_i - h*e_j)
              v[i] = vi - h;
              v[j] = vj - h;
              vsw = lambda_set.eval_value(v);

              // combine into difference quotient
              set_diff_quot(x, i, j, vne, vsw, vnw, vse, denom2);

              // restore coord
              v[j] = vj;
            }

            // restore coord
            v[i] = vi;
          }
        }
      }; // class LambdaHessHelper<...> (Richardson extrapolation)

      /**
       * \brief Lambda Function Evaluator class template
       *
       * This class implements the Evaluator class template of the actual Lambda function classes, because all of its
       * varying components have been outsourced to other helper classes, so the actual Evaluator class template is
       * syntactically identical for all lambda function implementations.
       *
       * \tparam LambdaSet_
       * The lambda set class defined in the lambda function class
       *
       * \tparam Traits_
       * The function evaluation traits, as defined in the Analytic::Function::Evaluator interface
       *
       * \tparam have_grad
       * \c true, if the user supplied all first-order partial derivatives as lambda expressions, otherwise \c false
       *
       * \tparam have_hess
       * \c true, if the user supplied all second-order partial derivatives as lambda expressions, otherwise \c false
       *
       *\ author Peter Zajac
       */
      template<typename LambdaSet_, typename Traits_, bool have_grad, bool have_hess>
      class LambdaFunctionEvaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        static constexpr int domain_dim = Traits_::domain_dim;
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;
        typedef typename Traits_::GradientType GradientType;
        typedef typename Traits_::HessianType HessianType;

        /// a local copy of the lambda set of this function
        LambdaSet_ lambda_set;

        /// gradient evaluation/extrapolation helper object
        Intern::LambdaGradHelper<LambdaSet_, domain_dim, DataType, ValueType, GradientType, have_grad> grad_helper;

        /// hessian evaluation/extrapolation helper object
        Intern::LambdaHessHelper<LambdaSet_, domain_dim, DataType, ValueType, HessianType, have_hess> hess_helper;

      public:
        template<typename LambdaFunction_>
        explicit LambdaFunctionEvaluator(const LambdaFunction_& function) :
          lambda_set(function.lambda_set),
          grad_helper(lambda_set, function.initial_h_grad, function.max_steps_grad),
          hess_helper(lambda_set, function.initial_h_hess, function.max_steps_hess)
        {
        }

        /// computes the function value
        ValueType value(const PointType& point)
        {
          return lambda_set.eval_value(point);
        }

        /// computes the function gradient
        GradientType gradient(const PointType& point)
        {
          GradientType grad;
          grad_helper.eval_grad(grad, point);
          return grad;
        }

        /// computes the function hessian
        HessianType hessian(const PointType& point)
        {
          HessianType hess;
          hess_helper.eval_hess(hess, point);
          return hess;
        }
      }; // class LambdaFunctionEvaluator<...>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Analytic 1D scalar lambda expression function implementation
     *
     * This class template acts as a wrapper class that implements the Analytic::Function interface for a
     * 1D scalar function by using C++ lambda expressions for the actual formulae evaluation.
     *
     * \see Please refer to the \ref analytic_lambda_function for details about how to use this class.
     *
     * \author Peter Zajac
     */
    template<typename LambdaValue_, typename LambdaDx_ = void, typename LambdaDxx_ = void>
    class LambdaScalarFunction1D :
      public Analytic::Function
    {
    public:
      /// this is a 1D function
      static constexpr int domain_dim = 1;
      /// this is a a scalar function
      typedef Analytic::Image::Scalar ImageType;

      /// we can compute values if LambdaValue_ is not void
      static constexpr bool have_value =  Intern::LambdaHelper<LambdaValue_>::have;
      /// we can compute gradients if LambdaDx_ and LambdaDy_ are not void
      static constexpr bool have_grad = Intern::LambdaHelper<LambdaDx_>::have;
      /// we can compute hessians if LambdaDxx_, LambdaDyy_ and LambdaDxy_ are not void
      static constexpr bool have_hess = Intern::LambdaHelper<LambdaDxx_>::have;

      /// nothing makes sense if we cannot compute function values
      static_assert(have_value, "LambdaScalarFunction1D needs at least the function value formula");

      /// we can always compute values
      static constexpr bool can_value =  true;
      /// we can always compute gradients (either directly or by Richardson extrapolation)
      static constexpr bool can_grad = true;
      /// we can always compute hessians (either directly or by Richardson extrapolation)
      static constexpr bool can_hess = true;

      /**
       * \brief Lambda function set container class
       *
       * This helper class is responsible for encapsulating all of the individual lambda functions helper classes
       * as well as offering a corresponding function value evaluation function which is used by the Richardson
       * extrapolation helper classes for gradient/hessian approximation.
       */
      class LambdaSet
      {
      public:
        /// the lambda for the function value
        Intern::LambdaHelper<LambdaValue_> value;
        /// the lambda for the X-derivative
        Intern::LambdaHelper<LambdaDx_> dx;
        /// the lambda for the XX-derivative
        Intern::LambdaHelper<LambdaDxx_> dxx;
        /// dummy helper for Y-derivative
        const Intern::LambdaHelper<void> dy;
        /// dummy helper for Z-derivative
        const Intern::LambdaHelper<void> dz;
        /// dummy helper for YY-derivative
        const Intern::LambdaHelper<void> dyy;
        /// dummy helper for ZZ-derivative
        const Intern::LambdaHelper<void> dzz;
        /// dummy helper for XY-derivative
        const Intern::LambdaHelper<void> dxy;
        /// dummy helper for YZ-derivative
        const Intern::LambdaHelper<void> dyz;
        /// dummy helper for ZX-derivative
        const Intern::LambdaHelper<void> dzx;

        /// constructor for values only
        explicit LambdaSet(LambdaValue_&& value_) :
          value(std::forward<LambdaValue_>(value_))
        {
        }

        /// constructor for values and gradients
        template<typename LDx_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_))
        {
        }

        /// constructor for values, gradients and hessians
        template<typename LDx_, typename LDxx_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_, LDxx_&& dxx_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_)),
          dxx(std::forward<LDxx_>(dxx_))
        {
        }

        /// function value evaluation function
        template<typename DT_>
        DT_ eval_value(const Tiny::Vector<DT_, 1>& v) const
        {
          return value.f(v[0]);
        }
      }; // class LambdaSet

      /// the actual lambda set for this function
      LambdaSet lambda_set;

      /// evaluator class template is outsourced
      template<typename Traits_>
      using Evaluator = Intern::LambdaFunctionEvaluator<LambdaSet, Traits_, have_grad, have_hess>;

      /// initial H for gradient Richardson extrapolation
      Real initial_h_grad = 0.01;
      /// initial H for hessian Richardson extrapolation
      Real initial_h_hess = 0.01;
      /// maximum number of gradient Richardson extrapolation steps
      int max_steps_grad = 10;
      /// maximum number of hessian Richardson extrapolation steps
      int max_steps_hess = 10;

      /**
       * \brief Constructor
       *
       * It is \b highly recommended to use the Analytic::create_lambda_function_scalar_1d() overloads to create objects of
       * this class. These overloads also offer a detailed description of the available parameters.
       */
      template<typename... Lambdas_>
      explicit LambdaScalarFunction1D(Lambdas_&&... lambdas) :
        lambda_set(std::forward<Lambdas_>(lambdas)...)
      {
      }
    }; // LambdaScalarFunction1D

    /**
     * \brief Creates a scalar 1D lambda function from function values only
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \returns
     * An instance of the LambdaScalarFunction1D class template
     */
    template<typename LambdaValue_>
    LambdaScalarFunction1D<LambdaValue_> create_lambda_function_scalar_1d(LambdaValue_&& value)
    {
      return LambdaScalarFunction1D<LambdaValue_>(std::forward<LambdaValue_>(value));
    }

    /**
     * \brief Creates a scalar 1D lambda function from function values and gradients
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx
     * The lambda expression for the first-order derivative of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction1D class template
     */
    template<typename LambdaValue_, typename LambdaDx_>
    LambdaScalarFunction1D<LambdaValue_, LambdaDx_> create_lambda_function_scalar_1d(LambdaValue_&& value, LambdaDx_&& dx)
    {
      return LambdaScalarFunction1D<LambdaValue_, LambdaDx_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx));
    }

    /**
     * \brief Creates a scalar 1D lambda function from function values, gradients and hessians
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx
     * The lambda expression for the first-order derivative of the scalar function.
     *
     * \param[in] dxx
     * The lambda expression for the second-order derivative of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction1D class template
     */
    template<typename LambdaValue_, typename LambdaDx_, typename LambdaDxx_>
    LambdaScalarFunction1D<LambdaValue_, LambdaDx_, LambdaDxx_> create_lambda_function_scalar_1d(
      LambdaValue_&& value, LambdaDx_&& dx, LambdaDxx_&& dxx)
    {
      return LambdaScalarFunction1D<LambdaValue_, LambdaDx_, LambdaDxx_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx), std::forward<LambdaDxx_>(dxx));
    }

    /**
     * \brief Analytic scalar 2D lambda expression function implementation
     *
     * This class template acts as a wrapper class that implements the Analytic::Function interface for a
     * 2D scalar function by using C++ lambda expressions for the actual formulae evaluation.
     *
     * \see Please refer to the \ref analytic_lambda_function for details about how to use this class.
     *
     * \author Peter Zajac
     */
    template<typename LambdaValue_,
      typename LambdaDx_ = void, typename LambdaDy_ = void,
      typename LambdaDxx_ = void, typename LambdaDyy_ = void, typename LambdaDxy_ = void>
    class LambdaScalarFunction2D :
      public Analytic::Function
    {
    public:
      /// this is a 2D function
      static constexpr int domain_dim = 2;
      /// this is a a scalar function
      typedef Analytic::Image::Scalar ImageType;

      /// we can compute values if LambdaValue_ is not void
      static constexpr bool have_value =  Intern::LambdaHelper<LambdaValue_>::have;
      /// we can compute gradients if LambdaDx_ and LambdaDy_ are not void
      static constexpr bool have_grad =
        Intern::LambdaHelper<LambdaDx_>::have &&
        Intern::LambdaHelper<LambdaDy_>::have;
      /// we can compute hessians if LambdaDxx_, LambdaDyy_ and LambdaDxy_ are not void
      static constexpr bool have_hess =
        Intern::LambdaHelper<LambdaDxx_>::have &&
        Intern::LambdaHelper<LambdaDyy_>::have &&
        Intern::LambdaHelper<LambdaDxy_>::have;

      /// nothing makes sense if we cannot compute function values
      static_assert(have_value, "LambdaScalarFunction2D needs at least the function value formula");

      /// we can always compute values
      static constexpr bool can_value =  true;
      /// we can always compute gradients (either directly or by Richardson extrapolation)
      static constexpr bool can_grad = true;
      /// we can always compute hessians (either directly or by Richardson extrapolation)
      static constexpr bool can_hess = true;

      /**
       * \brief Lambda function set container class
       *
       * This helper class is responsible for encapsulating all of the individual lambda functions helper classes
       * as well as offering a corresponding function value evaluation function which is used by the Richardson
       * extrapolation helper classes for gradient/hessian approximation.
       */
      class LambdaSet
      {
      public:
        /// the lambda for the function value
        Intern::LambdaHelper<LambdaValue_> value;
        /// the lambda for the X-derivative
        Intern::LambdaHelper<LambdaDx_> dx;
        /// the lambda for the Y-derivative
        Intern::LambdaHelper<LambdaDy_> dy;
        /// the lambda for the XX-derivative
        Intern::LambdaHelper<LambdaDxx_> dxx;
        /// the lambda for the YY-derivative
        Intern::LambdaHelper<LambdaDyy_> dyy;
        /// the lambda for the XY-derivative
        Intern::LambdaHelper<LambdaDxy_> dxy;
        /// dummy helper for Z-derivative
        const Intern::LambdaHelper<void> dz;
        /// dummy helper for ZZ-derivative
        const Intern::LambdaHelper<void> dzz;
        /// dummy helper for YZ-derivative
        const Intern::LambdaHelper<void> dyz;
        /// dummy helper for ZX-derivative
        const Intern::LambdaHelper<void> dzx;

        /// constructor for values only
        explicit LambdaSet(LambdaValue_&& value_) :
          value(std::forward<LambdaValue_>(value_))
        {
        }

        /// constructor for values and gradients
        template<typename LDx_, typename LDy_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_, LDy_&& dy_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_)),
          dy(std::forward<LDy_>(dy_))
        {
        }

        /// constructor for values, gradients and hessians
        template<typename LDx_, typename LDy_, typename LDxx_, typename LDyy_, typename LDxy_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_, LDy_&& dy_, LDxx_&& dxx_, LDyy_&& dyy_, LDxy_&& dxy_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_)),
          dy(std::forward<LDy_>(dy_)),
          dxx(std::forward<LDxx_>(dxx_)),
          dyy(std::forward<LDyy_>(dyy_)),
          dxy(std::forward<LDxy_>(dxy_))
        {
        }

        /// function value evaluation function
        template<typename DT_>
        DT_ eval_value(const Tiny::Vector<DT_, 2>& v) const
        {
          return value.f(v[0], v[1]);
        }
      }; // class LambdaSet

      /// the actual lambda set for this function
      LambdaSet lambda_set;

      /// evaluator class template is outsourced
      template<typename Traits_>
      using Evaluator = Intern::LambdaFunctionEvaluator<LambdaSet, Traits_, have_grad, have_hess>;

      /// initial H for gradient Richardson extrapolation
      Real initial_h_grad = 0.01;
      /// initial H for hessian Richardson extrapolation
      Real initial_h_hess = 0.01;
      /// maximum number of gradient Richardson extrapolation steps
      int max_steps_grad = 10;
      /// maximum number of hessian Richardson extrapolation steps
      int max_steps_hess = 10;

      /**
       * \brief Constructor
       *
       * It is \b highly recommended to use the Analytic::create_lambda_function_scalar_3d() overloads to create objects of
       * this class. These overloads also offer a detailed description of the available parameters.
       */
      template<typename... Lambdas_>
      explicit LambdaScalarFunction2D(Lambdas_&&... lambdas) :
        lambda_set(std::forward<Lambdas_>(lambdas)...)
      {
      }
    }; // LambdaScalarFunction2D

    /**
     * \brief Creates a scalar 2D lambda function from function values only
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \returns
     * An instance of the LambdaScalarFunction2D class template
     */
    template<typename LambdaValue_>
    LambdaScalarFunction2D<LambdaValue_> create_lambda_function_scalar_2d(LambdaValue_&& value)
    {
      return LambdaScalarFunction2D<LambdaValue_>(std::forward<LambdaValue_>(value));
    }

    /**
     * \brief Creates a scalar 2D lambda function from function values and gradients
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx, dy
     * The lambda expressions for the first-order partial derivatives of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction2D class template
     */
    template<typename LambdaValue_, typename LambdaDx_, typename LambdaDy_>
    LambdaScalarFunction2D<LambdaValue_, LambdaDx_, LambdaDy_> create_lambda_function_scalar_2d(
      LambdaValue_&& value, LambdaDx_&& dx, LambdaDy_&& dy)
    {
      return LambdaScalarFunction2D<LambdaValue_, LambdaDx_, LambdaDy_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx), std::forward<LambdaDy_>(dy));
    }

    /**
     * \brief Creates a scalar 2D lambda function from function values, gradients and hessians
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx, dy
     * The lambda expressions for the first-order partial derivatives of the scalar function.
     *
     * \param[in] dxx, dyy, dxy
     * The lambda expressions for the second-order partial derivatives of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction2D class template
     */
    template<typename LambdaValue_, typename LambdaDx_, typename LambdaDy_,
      typename LambdaDxx_, typename LambdaDyy_, typename LambdaDxy_>
    LambdaScalarFunction2D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDxx_, LambdaDyy_, LambdaDxy_> create_lambda_function_scalar_2d(
      LambdaValue_&& value, LambdaDx_&& dx, LambdaDy_&& dy, LambdaDxx_&& dxx, LambdaDyy_&& dyy, LambdaDxy_&& dxy)
    {
      return LambdaScalarFunction2D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDxx_, LambdaDyy_, LambdaDxy_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx), std::forward<LambdaDy_>(dy),
        std::forward<LambdaDxx_>(dxx), std::forward<LambdaDyy_>(dyy),std::forward<LambdaDxy_>(dxy));
    }

    /**
     * \brief Analytic 3D scalar lambda expression function implementation
     *
     * This class template acts as a wrapper class that implements the Analytic::Function interface for a
     * 3D scalar function by using C++ lambda expressions for the actual formulae evaluation.
     *
     * \see Please refer to the \ref analytic_lambda_function for details about how to use this class.
     *
     * \author Peter Zajac
     */
    template<typename LambdaValue_,
      typename LambdaDx_ = void, typename LambdaDy_ = void, typename LambdaDz_ = void,
      typename LambdaDxx_ = void, typename LambdaDyy_ = void, typename LambdaDzz_ = void,
      typename LambdaDxy_ = void, typename LambdaDyz_ = void, typename LambdaDzx_ = void>
    class LambdaScalarFunction3D :
      public Analytic::Function
    {
    public:
      /// this is a 3D function
      static constexpr int domain_dim = 3;
      /// this is a a scalar function
      typedef Analytic::Image::Scalar ImageType;

      /// we can compute values if LambdaValue_ is not void
      static constexpr bool have_value =  Intern::LambdaHelper<LambdaValue_>::have;
      /// are all first order lambdas given?
      static constexpr bool have_grad =
        Intern::LambdaHelper<LambdaDx_>::have &&
        Intern::LambdaHelper<LambdaDy_>::have &&
        Intern::LambdaHelper<LambdaDz_>::have;
      /// are all second order lambdas given?
      static constexpr bool have_hess =
        Intern::LambdaHelper<LambdaDxx_>::have &&
        Intern::LambdaHelper<LambdaDyy_>::have &&
        Intern::LambdaHelper<LambdaDzz_>::have &&
        Intern::LambdaHelper<LambdaDxy_>::have &&
        Intern::LambdaHelper<LambdaDyz_>::have &&
        Intern::LambdaHelper<LambdaDzx_>::have;

      /// nothing makes sense if we cannot compute function values
      static_assert(have_value, "LambdaScalarFunction3D needs at least the function value formula");

      /// we can always compute values
      static constexpr bool can_value =  true;
      /// we can always compute gradients (either directly or by Richardson extrapolation)
      static constexpr bool can_grad = true;
      /// we can always compute hessians (either directly or by Richardson extrapolation)
      static constexpr bool can_hess = true;

      /**
       * \brief Lambda function set container class
       *
       * This helper class is responsible for encapsulating all of the individual lambda functions helper classes
       * as well as offering a corresponding function value evaluation function which is used by the Richardson
       * extrapolation helper classes for gradient/hessian approximation.
       */
      class LambdaSet
      {
      public:
        /// the lambda for the function value
        Intern::LambdaHelper<LambdaValue_> value;
        /// the lambda for the X-derivative
        Intern::LambdaHelper<LambdaDx_> dx;
        /// the lambda for the Y-derivative
        Intern::LambdaHelper<LambdaDy_> dy;
        /// the lambda for the Z-derivative
        Intern::LambdaHelper<LambdaDz_> dz;
        /// the lambda for the XX-derivative
        Intern::LambdaHelper<LambdaDxx_> dxx;
        /// the lambda for the YY-derivative
        Intern::LambdaHelper<LambdaDyy_> dyy;
        /// the lambda for the ZZ-derivative
        Intern::LambdaHelper<LambdaDzz_> dzz;
        /// the lambda for the XY-derivative
        Intern::LambdaHelper<LambdaDxy_> dxy;
        /// the lambda for the YZ-derivative
        Intern::LambdaHelper<LambdaDyz_> dyz;
        /// the lambda for the ZX-derivative
        Intern::LambdaHelper<LambdaDzx_> dzx;

        /// constructor for values only
        explicit LambdaSet(LambdaValue_&& value_) :
          value(std::forward<LambdaValue_>(value_))
        {
        }

        /// constructor for values and gradients
        template<typename LDx_, typename LDy_, typename LDz_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_, LDy_&& dy_, LDz_&& dz_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_)),
          dy(std::forward<LDy_>(dy_)),
          dz(std::forward<LDz_>(dz_))
        {
        }

        /// constructor for values, gradients and hessians
        template<typename LDx_, typename LDy_, typename LDz_,
          typename LDxx_, typename LDyy_, typename LDzz_, typename LDxy_, typename LDyz_, typename LDzx_>
        explicit LambdaSet(LambdaValue_&& value_, LDx_&& dx_, LDy_&& dy_, LDz_&& dz_,
          LDxx_&& dxx_, LDyy_&& dyy_, LDzz_&& dzz_, LDxy_&& dxy_, LDyz_&& dyz_, LDzx_&& dzx_) :
          value(std::forward<LambdaValue_>(value_)),
          dx(std::forward<LDx_>(dx_)),
          dy(std::forward<LDy_>(dy_)),
          dz(std::forward<LDz_>(dz_)),
          dxx(std::forward<LDxx_>(dxx_)),
          dyy(std::forward<LDyy_>(dyy_)),
          dzz(std::forward<LDzz_>(dzz_)),
          dxy(std::forward<LDxy_>(dxy_)),
          dyz(std::forward<LDyz_>(dyz_)),
          dzx(std::forward<LDzx_>(dzx_))
        {
        }

        /// function value evaluation function
        template<typename DT_>
        DT_ eval_value(const Tiny::Vector<DT_, 3>& v) const
        {
          return value.f(v[0], v[1], v[2]);
        }
      }; // class LambdaSet

      /// the actual lambda set for this function
      LambdaSet lambda_set;

      /// evaluator class template is outsourced
      template<typename Traits_>
      using Evaluator = Intern::LambdaFunctionEvaluator<LambdaSet, Traits_, have_grad, have_hess>;

      /// initial H for gradient Richardson extrapolation
      Real initial_h_grad = 0.01;
      /// initial H for hessian Richardson extrapolation
      Real initial_h_hess = 0.01;
      /// maximum number of gradient Richardson extrapolation steps
      int max_steps_grad = 10;
      /// maximum number of hessian Richardson extrapolation steps
      int max_steps_hess = 10;

      /**
       * \brief Constructor
       *
       * It is \b highly recommended to use the Analytic::create_lambda_function_scalar_3d() overloads to create objects of
       * this class. These overloads also offer a detailed description of the available parameters.
       */
      template<typename... Lambdas_>
      explicit LambdaScalarFunction3D(Lambdas_&&... lambdas) :
        lambda_set(std::forward<Lambdas_>(lambdas)...)
      {
      }
    }; // LambdaScalarFunction3D

    /**
     * \brief Creates a scalar 3D lambda function from function values only
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \returns
     * An instance of the LambdaScalarFunction3D class template
     */
    template<typename LambdaValue_>
    LambdaScalarFunction3D<LambdaValue_> create_lambda_function_scalar_3d(LambdaValue_&& value)
    {
      return LambdaScalarFunction3D<LambdaValue_>(std::forward<LambdaValue_>(value));
    }

    /**
     * \brief Creates a scalar 3D lambda function from function values and gradients
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx, dy, dz
     * The lambda expressions for the first-order partial derivatives of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction3D class template
     */
    template<typename LambdaValue_, typename LambdaDx_, typename LambdaDy_, typename LambdaDz_>
    LambdaScalarFunction3D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDz_> create_lambda_function_scalar_3d(
      LambdaValue_&& value, LambdaDx_&& dx, LambdaDy_&& dy, LambdaDz_&& dz)
    {
      return LambdaScalarFunction3D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDz_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx), std::forward<LambdaDy_>(dy), std::forward<LambdaDz_>(dz));
    }

    /**
     * \brief Creates a 3D scalar lambda function from function values, gradients and hessians
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value
     * The lambda expression for the scalar function value.
     *
     * \param[in] dx, dy, dz
     * The lambda expressions for the first-order partial derivatives of the scalar function.
     *
     * \param[in] dxx, dyy, dzz, dxy, dyz, dzx
     * The lambda expressions for the second-order partial derivatives of the scalar function.
     *
     * \returns
     * An instance of the LambdaScalarFunction3D class template
     */
    template<typename LambdaValue_, typename LambdaDx_, typename LambdaDy_, typename LambdaDz_,
      typename LambdaDxx_, typename LambdaDyy_, typename LambdaDzz_,
      typename LambdaDxy_, typename LambdaDyz_, typename LambdaDzx_>
    LambdaScalarFunction3D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDz_,
      LambdaDxx_, LambdaDyy_, LambdaDzz_, LambdaDxy_, LambdaDyz_, LambdaDzx_> create_lambda_function_scalar_3d(
      LambdaValue_&& value, LambdaDx_&& dx, LambdaDy_&& dy, LambdaDz_&& dz,
      LambdaDxx_&& dxx, LambdaDyy_&& dyy, LambdaDzz_&& dzz, LambdaDxy_&& dxy, LambdaDyz_&& dyz, LambdaDzx_&& dzx)
    {
      return LambdaScalarFunction3D<LambdaValue_, LambdaDx_, LambdaDy_, LambdaDz_,
        LambdaDxx_, LambdaDyy_, LambdaDzz_, LambdaDxy_, LambdaDyz_, LambdaDzx_>(
        std::forward<LambdaValue_>(value), std::forward<LambdaDx_>(dx), std::forward<LambdaDy_>(dy), std::forward<LambdaDz_>(dz),
        std::forward<LambdaDxx_>(dxx), std::forward<LambdaDyy_>(dyy), std::forward<LambdaDzz_>(dzz),
        std::forward<LambdaDxy_>(dxy), std::forward<LambdaDyz_>(dyz), std::forward<LambdaDzx_>(dzx));
    }

    /**
     * \brief Analytic 2D vector-valued lambda expression function implementation
     *
     * This class template acts as a wrapper class that implements the Analytic::Function interface for a
     * 2D vector-valued function by using C++ lambda expressions for the actual formulae evaluation.
     *
     * \see Please refer to the \ref analytic_lambda_function for details about how to use this class.
     *
     * \author Peter Zajac
     */
    template<typename LambdaValue1_, typename LambdaValue2_,
      typename LambdaDx1_ = void, typename LambdaDx2_ = void,
      typename LambdaDy1_ = void, typename LambdaDy2_ = void,
      typename LambdaDxx1_ = void, typename LambdaDxx2_ = void,
      typename LambdaDyy1_ = void, typename LambdaDyy2_ = void,
      typename LambdaDxy1_ = void, typename LambdaDxy2_ = void>
    class LambdaVectorFunction2D :
      public Analytic::Function
    {
    public:
      /// this is a 2D function
      static constexpr int domain_dim = 2;
      /// this is a vector field
      typedef Analytic::Image::Vector<domain_dim> ImageType;

      /// we can compute values if LambdaValue_ is not void
      static constexpr bool have_value =
        Intern::LambdaHelper<LambdaValue1_>::have &&
        Intern::LambdaHelper<LambdaValue2_>::have;
      /// are all first order lambdas given?
      static constexpr bool have_grad =
        Intern::LambdaHelper<LambdaDx1_>::have &&
        Intern::LambdaHelper<LambdaDx2_>::have &&
        Intern::LambdaHelper<LambdaDy1_>::have &&
        Intern::LambdaHelper<LambdaDy2_>::have;
      /// are all second order lambdas given?
      static constexpr bool have_hess =
        Intern::LambdaHelper<LambdaDxx1_>::have &&
        Intern::LambdaHelper<LambdaDxx2_>::have &&
        Intern::LambdaHelper<LambdaDyy1_>::have &&
        Intern::LambdaHelper<LambdaDyy2_>::have &&
        Intern::LambdaHelper<LambdaDxy1_>::have &&
        Intern::LambdaHelper<LambdaDxy2_>::have;

      /// nothing makes sense if we cannot compute function values
      static_assert(have_value, "LambdaVectorFunction2D needs at least the function value formula");

      /// we can always compute values
      static constexpr bool can_value =  true;
      /// we can always compute gradients (either directly or by Richardson extrapolation)
      static constexpr bool can_grad = true;
      /// we can always compute hessians (either directly or by Richardson extrapolation)
      static constexpr bool can_hess = true;

      /**
       * \brief Lambda function set container class
       *
       * This helper class is responsible for encapsulating all of the individual lambda functions helper classes
       * as well as offering a corresponding function value evaluation function which is used by the Richardson
       * extrapolation helper classes for gradient/hessian approximation.
       */
      class LambdaSet
      {
      public:
        /// the lambdas for the function values
        Intern::LambdaHelper<LambdaValue1_> value1;
        Intern::LambdaHelper<LambdaValue2_> value2;
        Intern::LambdaHelper<void> value3;
        /// the lambdas for the X-derivatives
        Intern::LambdaHelper<LambdaDx1_> dx1;
        Intern::LambdaHelper<LambdaDx2_> dx2;
        Intern::LambdaHelper<void> dx3;
        /// the lambdas for the Y-derivatives
        Intern::LambdaHelper<LambdaDy1_> dy1;
        Intern::LambdaHelper<LambdaDy2_> dy2;
        Intern::LambdaHelper<void> dy3;
        /// the lambdas for the Z-derivatives
        Intern::LambdaHelper<void> dz1;
        Intern::LambdaHelper<void> dz2;
        Intern::LambdaHelper<void> dz3;
        /// the lambdas for the XX-derivatives
        Intern::LambdaHelper<LambdaDxx1_> dxx1;
        Intern::LambdaHelper<LambdaDxx2_> dxx2;
        Intern::LambdaHelper<void> dxx3;
        /// the lambdas for the YY-derivatives
        Intern::LambdaHelper<LambdaDyy1_> dyy1;
        Intern::LambdaHelper<LambdaDyy2_> dyy2;
        Intern::LambdaHelper<void> dyy3;
        /// the lambdas for the ZZ-derivatives
        Intern::LambdaHelper<void> dzz1;
        Intern::LambdaHelper<void> dzz2;
        Intern::LambdaHelper<void> dzz3;
        /// the lambdas for the XY-derivatives
        Intern::LambdaHelper<LambdaDxy1_> dxy1;
        Intern::LambdaHelper<LambdaDxy2_> dxy2;
        Intern::LambdaHelper<void> dxy3;
        /// the lambdas for the YZ-derivatives
        Intern::LambdaHelper<void> dyz1;
        Intern::LambdaHelper<void> dyz2;
        Intern::LambdaHelper<void> dyz3;
        /// the lambdas for the ZX-derivatives
        Intern::LambdaHelper<void> dzx1;
        Intern::LambdaHelper<void> dzx2;
        Intern::LambdaHelper<void> dzx3;

        /// constructor for values only
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_))
        {
        }

        /// constructor for values and gradients
        template<
          typename LDx1_, typename LDx2_,
          typename LDy1_, typename LDy2_>
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_,
          LDx1_&& dx1_, LDx2_&& dx2_, LDy1_&& dy1_, LDy2_&& dy2_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_)),
          dx1(std::forward<LDx1_>(dx1_)), dx2(std::forward<LDx2_>(dx2_)),
          dy1(std::forward<LDy1_>(dy1_)), dy2(std::forward<LDy2_>(dy2_))
        {
        }

        /// constructor for values, gradients and hessians
        template<
          typename LDx1_, typename LDx2_,
          typename LDy1_, typename LDy2_,
          typename LDxx1_, typename LDxx2_,
          typename LDyy1_, typename LDyy2_,
          typename LDxy1_, typename LDxy2_>
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_,
          LDx1_&& dx1_, LDx2_&& dx2_,
          LDy1_&& dy1_, LDy2_&& dy2_,
          LDxx1_&& dxx1_, LDxx2_&& dxx2_,
          LDyy1_&& dyy1_, LDyy2_&& dyy2_,
          LDxy1_&& dxy1_, LDxy2_&& dxy2_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_)),
          dx1(std::forward<LDx1_>(dx1_)), dx2(std::forward<LDx2_>(dx2_)),
          dy1(std::forward<LDy1_>(dy1_)), dy2(std::forward<LDy2_>(dy2_)),
          dxx1(std::forward<LDxx1_>(dxx1_)), dxx2(std::forward<LDxx2_>(dxx2_)),
          dyy1(std::forward<LDyy1_>(dyy1_)), dyy2(std::forward<LDyy2_>(dyy2_)),
          dxy1(std::forward<LDxy1_>(dxy1_)), dxy2(std::forward<LDxy2_>(dxy2_))
        {
        }

        /// function value evaluation function
        template<typename DT_>
        Tiny::Vector<DT_, 2> eval_value(const Tiny::Vector<DT_, 2>& v) const
        {
          return Tiny::Vector<DT_, 2>({value1.f(v[0], v[1]), value2.f(v[0], v[1])});
        }
      }; // class LambdaSet

      /// the actual lambda set for this function
      LambdaSet lambda_set;

      /// evaluator class template is outsourced
      template<typename Traits_>
      using Evaluator = Intern::LambdaFunctionEvaluator<LambdaSet, Traits_, have_grad, have_hess>;

      /// initial H for gradient Richardson extrapolation
      Real initial_h_grad = 0.01;
      /// initial H for hessian Richardson extrapolation
      Real initial_h_hess = 0.01;
      /// maximum number of gradient Richardson extrapolation steps
      int max_steps_grad = 10;
      /// maximum number of hessian Richardson extrapolation steps
      int max_steps_hess = 10;

      /**
       * \brief Constructor
       *
       * It is \b highly recommended to use the Analytic::create_lambda_function_vector_2d() overloads to create objects of
       * this class. These overloads also offer a detailed description of the available parameters.
       */
      template<typename... Lambdas_>
      explicit LambdaVectorFunction2D(Lambdas_&&... lambdas) :
        lambda_set(std::forward<Lambdas_>(lambdas)...)
      {
      }
    }; // LambdaVectorFunction2D

    /**
     * \brief Creates a vector-valued 2D lambda function from function values only
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction2D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_>
    LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_>
      create_lambda_function_vector_2d(LambdaValue1_&& value1, LambdaValue2_&& value2)
    {
      return LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_>(
        std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2));
    }

    /**
     * \brief Creates a vector-valued 2D lambda function from function values and gradients
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \param[in] dx1, dx2, dy1, dy2
     * The lambda expressions for the first-order partial derivatives of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction2D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_,
      typename LambdaDx1_, typename LambdaDx2_, typename LambdaDy1_, typename LambdaDy2_>
    LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_, LambdaDx1_, LambdaDx2_, LambdaDy1_, LambdaDy2_>
      create_lambda_function_vector_2d(LambdaValue1_&& value1, LambdaValue2_&& value2,
        LambdaDx1_&& dx1, LambdaDx2_&& dx2, LambdaDy1_&& dy1, LambdaDy2_&& dy2)
    {
      return LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_,
        LambdaDx1_, LambdaDx2_, LambdaDy1_, LambdaDy2_>(
          std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2),
          std::forward<LambdaDx1_>(dx1), std::forward<LambdaDx2_>(dx2),
          std::forward<LambdaDy1_>(dy1), std::forward<LambdaDy2_>(dy2));
    }

    /**
     * \brief Creates a vector-valued 2D lambda function from function values, gradients and hessians
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \param[in] dx1, dx2, dy1, dy2
     * The lambda expressions for the first-order partial derivatives of the vector field components.
     *
     * \param[in] dxx1, dxx2, dyy1, dyy2, dxy1, dxy2
     * The lambda expressions for the second-order partial derivatives of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction2D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_,
      typename LambdaDx1_, typename LambdaDx2_,
      typename LambdaDy1_, typename LambdaDy2_,
      typename LambdaDxx1_, typename LambdaDxx2_,
      typename LambdaDyy1_, typename LambdaDyy2_,
      typename LambdaDxy1_, typename LambdaDxy2_>
    LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_, LambdaDx1_, LambdaDx2_, LambdaDy1_, LambdaDy2_,
      LambdaDxx1_, LambdaDxx2_, LambdaDyy1_, LambdaDyy2_, LambdaDxy1_, LambdaDxy2_>
      create_lambda_function_vector_2d(LambdaValue1_&& value1, LambdaValue2_&& value2,
        LambdaDx1_&& dx1, LambdaDx2_&& dx2, LambdaDy1_&& dy1, LambdaDy2_&& dy2,
        LambdaDxx1_&& dxx1, LambdaDxx2_&& dxx2,
        LambdaDyy1_&& dyy1, LambdaDyy2_&& dyy2,
        LambdaDxy1_&& dxy1, LambdaDxy2_&& dxy2)
    {
      return LambdaVectorFunction2D<LambdaValue1_, LambdaValue2_, LambdaDx1_, LambdaDx2_, LambdaDy1_, LambdaDy2_,
        LambdaDxx1_, LambdaDxx2_, LambdaDyy1_, LambdaDyy2_, LambdaDxy1_, LambdaDxy2_>(
          std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2),
          std::forward<LambdaDx1_>(dx1), std::forward<LambdaDx2_>(dx2),
          std::forward<LambdaDy1_>(dy1), std::forward<LambdaDy2_>(dy2),
          std::forward<LambdaDxx1_>(dxx1), std::forward<LambdaDxx2_>(dxx2),
          std::forward<LambdaDyy1_>(dyy1), std::forward<LambdaDyy2_>(dyy2),
          std::forward<LambdaDxy1_>(dxy1), std::forward<LambdaDxy2_>(dxy2));
    }

    /**
     * \brief Analytic 3D vector-valued lambda expression function implementation
     *
     * This class template acts as a wrapper class that implements the Analytic::Function interface for a
     * 3D vector-valued function by using C++ lambda expressions for the actual formulae evaluation.
     *
     * \see Please refer to the \ref analytic_lambda_function for details about how to use this class.
     *
     * \author Peter Zajac
     */
    template<typename LambdaValue1_,typename LambdaValue2_,typename LambdaValue3_,
      typename LambdaDx1_ = void, typename LambdaDx2_ = void, typename LambdaDx3_ = void,
      typename LambdaDy1_ = void, typename LambdaDy2_ = void, typename LambdaDy3_ = void,
      typename LambdaDz1_ = void, typename LambdaDz2_ = void, typename LambdaDz3_ = void,
      typename LambdaDxx1_ = void, typename LambdaDxx2_ = void, typename LambdaDxx3_ = void,
      typename LambdaDyy1_ = void, typename LambdaDyy2_ = void, typename LambdaDyy3_ = void,
      typename LambdaDzz1_ = void, typename LambdaDzz2_ = void, typename LambdaDzz3_ = void,
      typename LambdaDxy1_ = void, typename LambdaDxy2_ = void, typename LambdaDxy3_ = void,
      typename LambdaDyz1_ = void, typename LambdaDyz2_ = void, typename LambdaDyz3_ = void,
      typename LambdaDzx1_ = void, typename LambdaDzx2_ = void, typename LambdaDzx3_ = void
    >
    class LambdaVectorFunction3D :
      public Analytic::Function
    {
    public:
      /// this is a 3D function
      static constexpr int domain_dim = 3;
      /// this is a vector field
      typedef Analytic::Image::Vector<domain_dim> ImageType;

      /// we can compute values if LambdaValue_ is not void
      static constexpr bool have_value =
        Intern::LambdaHelper<LambdaValue1_>::have &&
        Intern::LambdaHelper<LambdaValue2_>::have &&
        Intern::LambdaHelper<LambdaValue3_>::have;
        /// are all first order lambdas given?
      static constexpr bool have_grad =
        Intern::LambdaHelper<LambdaDx1_>::have &&
        Intern::LambdaHelper<LambdaDx2_>::have &&
        Intern::LambdaHelper<LambdaDx3_>::have &&
        Intern::LambdaHelper<LambdaDy1_>::have &&
        Intern::LambdaHelper<LambdaDy2_>::have &&
        Intern::LambdaHelper<LambdaDy3_>::have &&
        Intern::LambdaHelper<LambdaDz1_>::have &&
        Intern::LambdaHelper<LambdaDz2_>::have &&
        Intern::LambdaHelper<LambdaDz3_>::have;
      /// are all second order lambdas given?
      static constexpr bool have_hess =
        Intern::LambdaHelper<LambdaDxx1_>::have &&
        Intern::LambdaHelper<LambdaDxx2_>::have &&
        Intern::LambdaHelper<LambdaDxx3_>::have &&
        Intern::LambdaHelper<LambdaDyy1_>::have &&
        Intern::LambdaHelper<LambdaDyy2_>::have &&
        Intern::LambdaHelper<LambdaDyy3_>::have &&
        Intern::LambdaHelper<LambdaDzz1_>::have &&
        Intern::LambdaHelper<LambdaDzz2_>::have &&
        Intern::LambdaHelper<LambdaDzz3_>::have &&
        Intern::LambdaHelper<LambdaDxy1_>::have &&
        Intern::LambdaHelper<LambdaDxy2_>::have &&
        Intern::LambdaHelper<LambdaDxy3_>::have &&
        Intern::LambdaHelper<LambdaDyz1_>::have &&
        Intern::LambdaHelper<LambdaDyz2_>::have &&
        Intern::LambdaHelper<LambdaDyz3_>::have &&
        Intern::LambdaHelper<LambdaDzx1_>::have &&
        Intern::LambdaHelper<LambdaDzx2_>::have &&
        Intern::LambdaHelper<LambdaDzx3_>::have;

      /// nothing makes sense if we cannot compute function values
      static_assert(have_value, "LambdaVectorFunction3D needs at least the function value formula");

      /// we can always compute values
      static constexpr bool can_value =  true;
      /// we can always compute gradients (either directly or by Richardson extrapolation)
      static constexpr bool can_grad = true;
      /// we can always compute hessians (either directly or by Richardson extrapolation)
      static constexpr bool can_hess = true;

      /**
       * \brief Lambda function set container class
       *
       * This helper class is responsible for encapsulating all of the individual lambda functions helper classes
       * as well as offering a corresponding function value evaluation function which is used by the Richardson
       * extrapolation helper classes for gradient/hessian approximation.
       */
      class LambdaSet
      {
      public:
        /// the lambdas for the function values
        Intern::LambdaHelper<LambdaValue1_> value1;
        Intern::LambdaHelper<LambdaValue2_> value2;
        Intern::LambdaHelper<LambdaValue3_> value3;
        /// the lambdas for the X-derivatives
        Intern::LambdaHelper<LambdaDx1_> dx1;
        Intern::LambdaHelper<LambdaDx2_> dx2;
        Intern::LambdaHelper<LambdaDx3_> dx3;
        /// the lambdas for the Y-derivatives
        Intern::LambdaHelper<LambdaDy1_> dy1;
        Intern::LambdaHelper<LambdaDy2_> dy2;
        Intern::LambdaHelper<LambdaDy3_> dy3;
        /// the lambdas for the Z-derivatives
        Intern::LambdaHelper<LambdaDz1_> dz1;
        Intern::LambdaHelper<LambdaDz2_> dz2;
        Intern::LambdaHelper<LambdaDz3_> dz3;
        /// the lambdas for the XX-derivatives
        Intern::LambdaHelper<LambdaDxx1_> dxx1;
        Intern::LambdaHelper<LambdaDxx2_> dxx2;
        Intern::LambdaHelper<LambdaDxx3_> dxx3;
        /// the lambdas for the YY-derivatives
        Intern::LambdaHelper<LambdaDyy1_> dyy1;
        Intern::LambdaHelper<LambdaDyy2_> dyy2;
        Intern::LambdaHelper<LambdaDyy3_> dyy3;
        /// the lambdas for the ZZ-derivatives
        Intern::LambdaHelper<LambdaDzz1_> dzz1;
        Intern::LambdaHelper<LambdaDzz2_> dzz2;
        Intern::LambdaHelper<LambdaDzz3_> dzz3;
        /// the lambdas for the XY-derivatives
        Intern::LambdaHelper<LambdaDxy1_> dxy1;
        Intern::LambdaHelper<LambdaDxy2_> dxy2;
        Intern::LambdaHelper<LambdaDxy3_> dxy3;
        /// the lambdas for the YZ-derivatives
        Intern::LambdaHelper<LambdaDyz1_> dyz1;
        Intern::LambdaHelper<LambdaDyz2_> dyz2;
        Intern::LambdaHelper<LambdaDyz3_> dyz3;
        /// the lambdas for the ZX-derivatives
        Intern::LambdaHelper<LambdaDzx1_> dzx1;
        Intern::LambdaHelper<LambdaDzx2_> dzx2;
        Intern::LambdaHelper<LambdaDzx3_> dzx3;

        /// constructor for values only
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_, LambdaValue3_&& value3_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_)),
          value3(std::forward<LambdaValue3_>(value3_))
        {
        }

        /// constructor for values and gradients
        template<
          typename LDx1_, typename LDx2_, typename LDx3_,
          typename LDy1_, typename LDy2_, typename LDy3_,
          typename LDz1_, typename LDz2_, typename LDz3_>
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_, LambdaValue3_&& value3_,
          LDx1_&& dx1_, LDx2_&& dx2_, LDx3_&& dx3_,
          LDy1_&& dy1_, LDy2_&& dy2_, LDy3_&& dy3_,
          LDz1_&& dz1_, LDz2_&& dz2_, LDz3_&& dz3_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_)),
          value3(std::forward<LambdaValue3_>(value3_)),
          dx1(std::forward<LDx1_>(dx1_)), dx2(std::forward<LDx2_>(dx2_)), dx3(std::forward<LDx3_>(dx3_)),
          dy1(std::forward<LDy1_>(dy1_)), dy2(std::forward<LDy2_>(dy2_)), dy3(std::forward<LDy3_>(dy3_)),
          dz1(std::forward<LDz1_>(dz1_)), dz2(std::forward<LDz2_>(dz2_)), dz3(std::forward<LDz3_>(dz3_))
        {
        }

        /// constructor for values, gradients and hessians
        template<
          typename LDx1_, typename LDx2_, typename LDx3_,
          typename LDy1_, typename LDy2_, typename LDy3_,
          typename LDz1_, typename LDz2_, typename LDz3_,
          typename LDxx1_, typename LDxx2_, typename LDxx3_,
          typename LDyy1_, typename LDyy2_, typename LDyy3_,
          typename LDzz1_, typename LDzz2_, typename LDzz3_,
          typename LDxy1_, typename LDxy2_, typename LDxy3_,
          typename LDyz1_, typename LDyz2_, typename LDyz3_,
          typename LDzx1_, typename LDzx2_, typename LDzx3_>
        explicit LambdaSet(LambdaValue1_&& value1_, LambdaValue2_&& value2_, LambdaValue3_&& value3_,
          LDx1_&& dx1_, LDx2_&& dx2_, LDx3_&& dx3_,
          LDy1_&& dy1_, LDy2_&& dy2_, LDy3_&& dy3_,
          LDz1_&& dz1_, LDz2_&& dz2_, LDz3_&& dz3_,
          LDxx1_&& dxx1_, LDxx2_&& dxx2_, LDxx3_&& dxx3_,
          LDyy1_&& dyy1_, LDyy2_&& dyy2_, LDyy3_&& dyy3_,
          LDzz1_&& dzz1_, LDzz2_&& dzz2_, LDzz3_&& dzz3_,
          LDxy1_&& dxy1_, LDxy2_&& dxy2_, LDxy3_&& dxy3_,
          LDyz1_&& dyz1_, LDyz2_&& dyz2_, LDyz3_&& dyz3_,
          LDzx1_&& dzx1_, LDzx2_&& dzx2_, LDzx3_&& dzx3_) :
          value1(std::forward<LambdaValue1_>(value1_)),
          value2(std::forward<LambdaValue2_>(value2_)),
          value3(std::forward<LambdaValue3_>(value3_)),
          dx1(std::forward<LDx1_>(dx1_)), dx2(std::forward<LDx2_>(dx2_)), dx3(std::forward<LDx3_>(dx3_)),
          dy1(std::forward<LDy1_>(dy1_)), dy2(std::forward<LDy2_>(dy2_)), dy3(std::forward<LDy3_>(dy3_)),
          dz1(std::forward<LDz1_>(dz1_)), dz2(std::forward<LDz2_>(dz2_)), dz3(std::forward<LDz3_>(dz3_)),
          dxx1(std::forward<LDxx1_>(dxx1_)), dxx2(std::forward<LDxx2_>(dxx2_)), dxx3(std::forward<LDxx3_>(dxx3_)),
          dyy1(std::forward<LDyy1_>(dyy1_)), dyy2(std::forward<LDyy2_>(dyy2_)), dyy3(std::forward<LDyy3_>(dyy3_)),
          dzz1(std::forward<LDzz1_>(dzz1_)), dzz2(std::forward<LDzz2_>(dzz2_)), dzz3(std::forward<LDzz3_>(dzz3_)),
          dxy1(std::forward<LDxy1_>(dxy1_)), dxy2(std::forward<LDxy2_>(dxy2_)), dxy3(std::forward<LDxy3_>(dxy3_)),
          dyz1(std::forward<LDyz1_>(dyz1_)), dyz2(std::forward<LDyz2_>(dyz2_)), dyz3(std::forward<LDyz3_>(dyz3_)),
          dzx1(std::forward<LDzx1_>(dzx1_)), dzx2(std::forward<LDzx2_>(dzx2_)), dzx3(std::forward<LDzx3_>(dzx3_))
        {
        }

        /// function value evaluation function
        template<typename DT_>
        Tiny::Vector<DT_, 3> eval_value(const Tiny::Vector<DT_, 3>& v) const
        {
          return Tiny::Vector<DT_, 3>({
            value1.f(v[0], v[1], v[2]),
            value2.f(v[0], v[1], v[2]),
            value3.f(v[0], v[1], v[2])});
        }
      }; // class LambdaSet

      /// the actual lambda set for this function
      LambdaSet lambda_set;

      /// evaluator class template is outsourced
      template<typename Traits_>
      using Evaluator = Intern::LambdaFunctionEvaluator<LambdaSet, Traits_, have_grad, have_hess>;

      /// initial H for gradient Richardson extrapolation
      Real initial_h_grad = 0.01;
      /// initial H for hessian Richardson extrapolation
      Real initial_h_hess = 0.01;
      /// maximum number of gradient Richardson extrapolation steps
      int max_steps_grad = 10;
      /// maximum number of hessian Richardson extrapolation steps
      int max_steps_hess = 10;

      /**
       * \brief Constructor
       *
       * It is \b highly recommended to use the Analytic::create_lambda_function_vector_3d() overloads to create objects of
       * this class. These overloads also offer a detailed description of the available parameters.
       */
      template<typename... Lambdas_>
      explicit LambdaVectorFunction3D(Lambdas_&&... lambdas) :
        lambda_set(std::forward<Lambdas_>(lambdas)...)
      {
      }
    }; // LambdaVectorFunction3D

    /**
     * \brief Creates a vector-valued 3D lambda function from function values only
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction3D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_, typename LambdaValue3_>
    LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_>
      create_lambda_function_vector_3d(LambdaValue1_&& value1, LambdaValue2_&& value2, LambdaValue3_&& value3)
    {
      return LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_>(
        std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2), std::forward<LambdaValue3_>(value3));
    }

    /**
     * \brief Creates a vector-valued 3D lambda function from function values and gradients
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \param[in] dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3
     * The lambda expressions for the first-order partial derivatives of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction3D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_, typename LambdaValue3_,
      typename LambdaDx1_, typename LambdaDx2_, typename LambdaDx3_,
      typename LambdaDy1_, typename LambdaDy2_, typename LambdaDy3_,
      typename LambdaDz1_, typename LambdaDz2_, typename LambdaDz3_>
    LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_,
      LambdaDx1_, LambdaDx2_, LambdaDx3_, LambdaDy1_, LambdaDy2_, LambdaDy3_, LambdaDz1_, LambdaDz2_, LambdaDz3_>
      create_lambda_function_vector_3d(LambdaValue1_&& value1, LambdaValue2_&& value2, LambdaValue3_&& value3,
        LambdaDx1_&& dx1, LambdaDx2_&& dx2, LambdaDx3_&& dx3,
        LambdaDy1_&& dy1, LambdaDy2_&& dy2, LambdaDy3_&& dy3,
        LambdaDz1_&& dz1, LambdaDz2_&& dz2, LambdaDz3_&& dz3)
    {
      return LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_,
        LambdaDx1_, LambdaDx2_, LambdaDx3_, LambdaDy1_, LambdaDy2_, LambdaDy3_, LambdaDz1_, LambdaDz2_, LambdaDz3_>(
        std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2), std::forward<LambdaValue3_>(value3),
        std::forward<LambdaDx1_>(dx1), std::forward<LambdaDx2_>(dx2), std::forward<LambdaDx3_>(dx3),
        std::forward<LambdaDy1_>(dy1), std::forward<LambdaDy2_>(dy2), std::forward<LambdaDy3_>(dy3),
        std::forward<LambdaDz1_>(dz1), std::forward<LambdaDz2_>(dz2), std::forward<LambdaDz3_>(dz3));
    }

    /**
     * \brief Creates a vector-valued 3D lambda function from function values, gradients and hessians
     *
     * \see See \ref analytic_lambda_function page for details.
     *
     * \param[in] value1, value2, value3
     * The lambda expressions for the function values of the vector field components.
     *
     * \param[in] dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3
     * The lambda expressions for the first-order partial derivatives of the vector field components.
     *
     * \param[in] dxx1, dxx2, dxx3, dyy1, dyy2, dyy3, dzz1, dzz2, dzz3, dxy1, dxy2, dxy3, dyz1, dyz2, dyz3, dzx1, dzx2, dzx3
     * The lambda expressions for the second-order partial derivatives of the vector field components.
     *
     * \returns
     * An instance of the LambdaVectorFunction3D class template
     */
    template<typename LambdaValue1_, typename LambdaValue2_, typename LambdaValue3_,
      typename LambdaDx1_, typename LambdaDx2_, typename LambdaDx3_,
      typename LambdaDy1_, typename LambdaDy2_, typename LambdaDy3_,
      typename LambdaDz1_, typename LambdaDz2_, typename LambdaDz3_,
      typename LambdaDxx1_, typename LambdaDxx2_, typename LambdaDxx3_,
      typename LambdaDyy1_, typename LambdaDyy2_, typename LambdaDyy3_,
      typename LambdaDzz1_, typename LambdaDzz2_, typename LambdaDzz3_,
      typename LambdaDxy1_, typename LambdaDxy2_, typename LambdaDxy3_,
      typename LambdaDyz1_, typename LambdaDyz2_, typename LambdaDyz3_,
      typename LambdaDzx1_, typename LambdaDzx2_, typename LambdaDzx3_>
    LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_,
      LambdaDx1_, LambdaDx2_, LambdaDx3_, LambdaDy1_, LambdaDy2_, LambdaDy3_, LambdaDz1_, LambdaDz2_, LambdaDz3_,
      LambdaDxx1_, LambdaDxx2_, LambdaDxx3_, LambdaDyy1_, LambdaDyy2_, LambdaDyy3_, LambdaDzz1_, LambdaDzz2_, LambdaDzz3_,
      LambdaDxy1_, LambdaDxy2_, LambdaDxy3_, LambdaDyz1_, LambdaDyz2_, LambdaDyz3_, LambdaDzx1_, LambdaDzx2_, LambdaDzx3_>
      create_lambda_function_vector_3d(LambdaValue1_&& value1, LambdaValue2_&& value2, LambdaValue3_&& value3,
      LambdaDx1_&& dx1, LambdaDx2_&& dx2, LambdaDx3_&& dx3,
      LambdaDy1_&& dy1, LambdaDy2_&& dy2, LambdaDy3_&& dy3,
      LambdaDz1_&& dz1, LambdaDz2_&& dz2, LambdaDz3_&& dz3,
      LambdaDxx1_&& dxx1, LambdaDxx2_&& dxx2, LambdaDxx3_&& dxx3,
      LambdaDyy1_&& dyy1, LambdaDyy2_&& dyy2, LambdaDyy3_&& dyy3,
      LambdaDzz1_&& dzz1, LambdaDzz2_&& dzz2, LambdaDzz3_&& dzz3,
      LambdaDxy1_&& dxy1, LambdaDxy2_&& dxy2, LambdaDxy3_&& dxy3,
      LambdaDyz1_&& dyz1, LambdaDyz2_&& dyz2, LambdaDyz3_&& dyz3,
      LambdaDzx1_&& dzx1, LambdaDzx2_&& dzx2, LambdaDzx3_&& dzx3)
    {
      return LambdaVectorFunction3D<LambdaValue1_, LambdaValue2_, LambdaValue3_,
        LambdaDx1_, LambdaDx2_, LambdaDx3_, LambdaDy1_, LambdaDy2_, LambdaDy3_, LambdaDz1_, LambdaDz2_, LambdaDz3_,
        LambdaDxx1_, LambdaDxx2_, LambdaDxx3_, LambdaDyy1_, LambdaDyy2_, LambdaDyy3_, LambdaDzz1_, LambdaDzz2_, LambdaDzz3_,
        LambdaDxy1_, LambdaDxy2_, LambdaDxy3_, LambdaDyz1_, LambdaDyz2_, LambdaDyz3_, LambdaDzx1_, LambdaDzx2_, LambdaDzx3_>(
        std::forward<LambdaValue1_>(value1), std::forward<LambdaValue2_>(value2), std::forward<LambdaValue3_>(value3),
        std::forward<LambdaDx1_>(dx1), std::forward<LambdaDx2_>(dx2), std::forward<LambdaDx3_>(dx3),
        std::forward<LambdaDy1_>(dy1), std::forward<LambdaDy2_>(dy2), std::forward<LambdaDy3_>(dy3),
        std::forward<LambdaDz1_>(dz1), std::forward<LambdaDz2_>(dz2), std::forward<LambdaDz3_>(dz3),
        std::forward<LambdaDxx1_>(dxx1), std::forward<LambdaDxx2_>(dxx2), std::forward<LambdaDxx3_>(dxx3),
        std::forward<LambdaDyy1_>(dyy1), std::forward<LambdaDyy2_>(dyy2), std::forward<LambdaDyy3_>(dyy3),
        std::forward<LambdaDzz1_>(dzz1), std::forward<LambdaDzz2_>(dzz2), std::forward<LambdaDzz3_>(dzz3),
        std::forward<LambdaDxy1_>(dxy1), std::forward<LambdaDxy2_>(dxy2), std::forward<LambdaDxy3_>(dxy3),
        std::forward<LambdaDyz1_>(dyz1), std::forward<LambdaDyz2_>(dyz2), std::forward<LambdaDyz3_>(dyz3),
        std::forward<LambdaDzx1_>(dzx1), std::forward<LambdaDzx2_>(dzx2), std::forward<LambdaDzx3_>(dzx3));
    }
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_LAMBDA_FUNCTION_HPP
