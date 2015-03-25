#pragma once
#ifndef KERNEL_SPACE_HERMITE3_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_HERMITE3_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Hermite3
    {
      template<
        typename Space_,
        typename Function_,
        typename Shape_,
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
      class NodeFunctional<Space_, Function_, Shape::Simplex<2>, 0, DataType_> :
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
          static constexpr bool need_value = true;
          static constexpr bool need_grad = true;
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;
        typedef typename Function_::template Evaluator<AnalyticEvalTraits> FuncEval;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;

        typename AnalyticEvalTraits::ValueType _func_value;
        typename AnalyticEvalTraits::GradientType _func_grad;

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
        }

        void finish()
        {
          _func_eval.finish();
          _trafo_eval.finish();
          BaseClass::finish();
        }

        Index get_max_assigned_dofs() const
        {
          return 3;
        }

        Index get_num_assigned_dofs() const
        {
          return 3;
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
          default:
            return DataType_(0);
          }
        }
      };

      template<
        typename Space_,
        typename Function_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, Shape::Simplex<2>, 2, DataType_> :
        public NodeFunctionalBase<Space_, Function_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Function_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<Shape::Simplex<2>, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct FunctionConfig :
          public Trafo::AnalyticConfigBase
        {
          static constexpr bool need_value = true;
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

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
          dom_point = DataType_(1) / DataType_(3);

          // compute trafo-data
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          return _func_eval.value(trafo_data);
        }
      };

      template<
        typename Space_,
        typename Function_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, Shape::Hypercube<1>, 0, DataType_> :
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
          static constexpr bool need_value = true;
          static constexpr bool need_grad = true;
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;
        typedef typename Function_::template Evaluator<AnalyticEvalTraits> FuncEval;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;

        typename AnalyticEvalTraits::ValueType _func_value;
        typename AnalyticEvalTraits::GradientType _func_grad;

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
        }

        void finish()
        {
          _func_eval.finish();
          _trafo_eval.finish();
          BaseClass::finish();
        }

        Index get_max_assigned_dofs() const
        {
          return 2;
        }

        Index get_num_assigned_dofs() const
        {
          return 2;
        }

        DataType_ operator()(Index assign_idx) const
        {
          switch(int(assign_idx))
          {
          case 0:
            return _func_value;
          case 1:
            return _func_grad(Index(0));
          default:
            return DataType_(0);
          }
        }
      };

      template<
        typename Space_,
        typename Function_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, Shape::Hypercube<2>, 0, DataType_> :
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
          static constexpr bool need_value = true;
          static constexpr bool need_grad = true;
          static constexpr bool need_hess = true;
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
          return 4;
        }

        Index get_num_assigned_dofs() const
        {
          return 4;
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
            return _func_hess(Index(0), Index(1));
          default:
            return DataType_(0);
          }
        }
      };
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_HERMITE3_NODE_FUNCTIONAL_HPP
