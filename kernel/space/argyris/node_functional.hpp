#pragma once
#ifndef KERNEL_SPACE_ARGYRIS_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_ARGYRIS_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Argyris
    {
      template<
        typename Space_,
        typename Function_,
        int dim_,
        typename DataType_>
      class NodeFunctional :
        public NodeFunctionalNull<Space_, Function_, DataType_>
      {
      public:
        explicit NodeFunctional(const Space_& space, const Function_& function) :
          NodeFunctionalNull<Space_, Function_, DataType_>(space, function)
        {
        }
      };

      template<
        typename Space_,
        typename Function_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, 0, DataType_> :
        public NodeFunctionalBase<Space_, Function_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Function_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<Shape::Vertex, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct FunctionConfig :
          public Trafo::AnalyticConfigBase
        {
          enum
          {
            need_value = 1,
            need_grad = 1,
            need_hess = 1
          };
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;
        typedef typename Function_::template Evaluator<AnalyticEvalTraits> FuncEval;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;

        typename AnalyticEvalTraits::ValueType _func_value;
        typename AnalyticEvalTraits::GradientType _func_grad;
        typename AnalyticEvalTraits::HessianType _func_hess;

      public:
        explicit NodeFunctional(const Space_& space, const Function_& function) :
          BaseClass(space, function),
          _trafo_eval(space.get_trafo()),
          _func_eval(function)
        {
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);

          // prepare
          _trafo_eval.prepare(cell_index);
          _func_eval.prepare(_trafo_eval);

          // compute trafo-data
          DomainPointType dom_point;
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          _func_value = _func_eval.value(trafo_data);
          _func_grad = _func_eval.gradient(trafo_data);
          _func_hess = _func_eval.hessian(trafo_data);
        }

        void finish()
        {
          _func_eval.finish();
          _trafo_eval.finish();
          BaseClass::finish();
        }

        Index get_max_assigned_dofs() const
        {
          return 6;
        }

        Index get_num_assigned_dofs() const
        {
          return 6;
        }

        DataType_ operator()(Index assign_idx) const
        {
          switch(int(assign_idx))
          {
          case 0:
            return _func_value;
          case 1:
            return _func_grad(Index(0));
          case 2:
            return _func_grad(Index(1));
          case 3:
            return _func_hess(Index(0), Index(0));
          case 4:
            return _func_hess(Index(1), Index(1));
          case 5:
            return _func_hess(Index(0), Index(1));
          default:
            return DataType_(0);
          }
        }
      };

      template<
        typename Space_,
        typename Function_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, 1, DataType_> :
        public NodeFunctionalBase<Space_, Function_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Function_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<Shape::Simplex<1>, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct FunctionConfig :
          public Trafo::AnalyticConfigBase
        {
          enum
          {
            need_grad = 1
          };
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig FunctionTrafoConfig;

        struct TrafoConfig :
          public FunctionTrafoConfig
        {
          enum
          {
            need_jac_mat = 1
          };
        };

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;
        typedef typename Function_::template Evaluator<AnalyticEvalTraits> FuncEval;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;

      public:
        explicit NodeFunctional(const Space_& space, const Function_& function) :
          BaseClass(space, function),
          _trafo_eval(space.get_trafo()),
          _func_eval(function)
        {
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);
          _func_eval.prepare(_trafo_eval);
        }

        void finish()
        {
          _func_eval.finish();
          _trafo_eval.finish();
          BaseClass::finish();
        }

        Index get_max_assigned_dofs() const
        {
          return 1;
        }

        Index get_num_assigned_dofs() const
        {
          return 1;
        }

        DataType_ operator()(Index /*assign_idx*/) const
        {
          DomainPointType dom_point;
          dom_point(0) = DataType_(0.5);

          // compute trafo-data
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // compute edge normal
          DataType_ dnx = +trafo_data.jac_mat(1,0);
          DataType_ dny = -trafo_data.jac_mat(0,0);
          DataType_ dnl = DataType_(1) / Math::sqrt(dnx*dnx + dny*dny);
          dnx *= dnl;
          dny *= dnl;

          // evaluate function gradient
          typename AnalyticEvalTraits::GradientType func_grad = _func_eval.gradient(trafo_data);

          // return edge normal derivative
          return dnx * func_grad(0) + dny * func_grad(1);
        }
      };
    } // namespace Argyris
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_ARGYRIS_NODE_FUNCTIONAL_HPP
