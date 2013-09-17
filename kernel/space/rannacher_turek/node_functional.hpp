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
        typename Functor_,
        int codim_,
        typename Variant_,
        typename DataType_>
      class NodeFunctional :
        public NodeFunctionalNull<Space_, Functor_, DataType_>
      {
      public:
        explicit NodeFunctional(const Space_& space, const Functor_& functor) :
          NodeFunctionalNull<Space_, Functor_, DataType_>(space, functor)
        {
        }
      };

      template<
        typename Space_,
        typename Functor_,
        typename Variant_,
        typename DataType_>
      class NodeFunctional<Space_, Functor_, 1, Variant_, DataType_> :
        public NodeFunctionalBase<Space_, Functor_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Functor_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;
        typedef typename TrafoType::template Evaluator<FacetType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_img_point = 1,
            need_jac_det = 1
          };
        };

        typedef typename TrafoEvalType::template TrafoConfig<TrafoConfig> TrafoDataConfig;

        typedef Trafo::EvalData<TrafoEvalTraits, TrafoDataConfig> TrafoEvalData;

        struct FuncEvalTraits
        {
          enum
          {
            image_dim = TrafoEvalTraits::image_dim
          };

          typedef TrafoEvalType TrafoEvaluator;
          typedef TrafoEvalData TrafoData;
          typedef DataType_ DataType;
          typedef DataType_ ValueType;
        };

        typedef typename Functor_::template ValueEvaluator<FuncEvalTraits> FuncEval;

        typedef Cubature::Rule<FacetType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space, const Functor_& functor) :
          BaseClass(space, functor),
          _trafo_eval(space.get_trafo()),
          _func_eval(functor),
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
          DataType_ value(DataType_(0)), mean(DataType_(0)), v;
          TrafoEvalData trafo_data;
          const Index n(_cub_rule.get_num_points());
          for(Index i(0); i < n; ++i)
          {
            trafo_data(_trafo_eval, _cub_rule.get_point(i));
            _func_eval(v, trafo_data);
            value += _cub_rule.get_weight(i) * trafo_data.jac_det * v;
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }
          return value / mean;
        }
      };
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_NODE_FUNCTIONAL_HPP
