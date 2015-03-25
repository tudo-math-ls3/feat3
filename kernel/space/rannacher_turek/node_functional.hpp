#pragma once
#ifndef KERNEL_SPACE_RANNACHER_TUREK_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_RANNACHER_TUREK_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace RannacherTurek
    {
      template<
        typename Space_,
        typename Function_,
        int codim_,
        typename Variant_,
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
        typename Variant_,
        typename DataType_>
      class NodeFunctional<Space_, Function_, 1, Variant_, DataType_> :
        public NodeFunctionalBase<Space_, Function_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Function_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;
        typedef typename TrafoType::template Evaluator<FacetType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct FunctionConfig :
          public Trafo::AnalyticConfigBase
        {
          static constexpr bool need_value = true;
          static constexpr bool need_grad = false;
          static constexpr bool need_hess = false;
        };

        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig FunctionTrafoConfig;

        struct TrafoConfig :
          public FunctionTrafoConfig
        {
          static constexpr bool need_jac_det = true;
        };

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;
        typedef typename Function_::template Evaluator<AnalyticEvalTraits> FuncEval;

        typedef Cubature::Rule<FacetType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space, const Function_& function) :
          BaseClass(space, function),
          _trafo_eval(space.get_trafo()),
          _func_eval(function),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:2"))
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
          DataType_ value(DataType_(0)), mean(DataType_(0));
          TrafoEvalData trafo_data;
          const Index n(_cub_rule.get_num_points());
          for(Index i(0); i < n; ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += _cub_rule.get_weight(i) * trafo_data.jac_det * _func_eval.value(trafo_data);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }
          return value / mean;
        }
      };
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_NODE_FUNCTIONAL_HPP
