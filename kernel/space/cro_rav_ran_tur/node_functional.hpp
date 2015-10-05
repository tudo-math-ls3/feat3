#pragma once
#ifndef KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace CroRavRanTur
    {
      template<
        typename Space_,
        typename Shape_,
        int codim_,
        typename DataType_>
      class NodeFunctional :
        public NodeFunctionalNull<Space_, DataType_>
      {
      public:
        explicit NodeFunctional(const Space_& space) :
          NodeFunctionalNull<Space_, DataType_>(space)
        {
        }
      };

      /**
       * \brief Node Functional implementation for Crouzeix-Raviart
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        int shape_dim_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Simplex<shape_dim_>, 1, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 1;

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

        TrafoEvalType _trafo_eval;
        DomainPointType _barycentre;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
        {
          // initialize facet barycentre
          for(int i(0); i < ShapeType::dimension-1; ++i)
            _barycentre[i] = DataType_(1) / DataType_(shape_dim_);
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);
        }

        void finish()
        {
          _trafo_eval.finish();
          BaseClass::finish();
        }

        int get_num_assigned_dofs() const
        {
          return max_assigned_dofs;
        }

        template<typename NodeData_, typename Function_>
        void operator()(NodeData_& node_data, const Function_& function) const
        {
          // fetch the function's trafo config
          typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

          // declare trafo evaluation data
          typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

          // declare evaluation traits
          typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

          // prepare function evaluator
          func_eval.prepare(_trafo_eval);

          // compute trafo data
          DomainPointType dom_point;
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, _barycentre);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data);

          // finish function evaluator
          func_eval.finish();
        }
      };

      /**
       * \brief Node Functional implementation for Rannacher-Turek
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        int shape_dim_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<shape_dim_>, 1, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 1;

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

        template<typename FuncTrafoConfig_>
        struct TrafoConfig :
          public FuncTrafoConfig_
        {
          static constexpr bool need_jac_det = true;
        };

        typedef Cubature::Rule<FacetType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:2"))
        {
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);
        }

        void finish()
        {
          _trafo_eval.finish();
          BaseClass::finish();
        }

        int get_num_assigned_dofs() const
        {
          return max_assigned_dofs;
        }

        template<typename NodeData_, typename Function_>
        void operator()(NodeData_& node_data, const Function_& function) const
        {
          // fetch the function's trafo config
          typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig FunctionTrafoConfig;

          // declare trafo evaluation data
          typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig<FunctionTrafoConfig>>::EvalDataType TrafoEvalData;

          // declare evaluation traits
          typedef Trafo::AnalyticEvalTraits<TrafoEvalType, TrafoEvalData> AnalyticEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

          // prepare function evaluator
          func_eval.prepare(_trafo_eval);

          typename AnalyticEvalTraits::ValueType value(DataType_(0));
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;

          // integrate over facet
          for(int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += _cub_rule.get_weight(i) * trafo_data.jac_det * func_eval.value(trafo_data);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set integral mean
          node_data[0] = (DataType_(1) / mean) * value;

          // finish function evaluator
          func_eval.finish();
        }
      };
    } // namespace CroRavRanTur
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP
