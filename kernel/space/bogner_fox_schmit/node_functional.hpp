#pragma once
#ifndef KERNEL_SPACE_BOGNER_FOX_SCHMIT_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_BOGNER_FOX_SCHMIT_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace BognerFoxSchmit
    {
      template<
        typename Space_,
        typename Functor_,
        int space_dim_,
        int shape_dim_,
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
        typename DataType_>
      class NodeFunctional<Space_, Functor_, 1, 0, DataType_> :
        public NodeFunctionalBase<Space_, Functor_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Functor_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename TrafoType::template Evaluator<Shape::Vertex, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_img_point = 1
          };
        };

        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        struct FuncValueEvalTraits
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

        typedef typename Functor_::template ValueEvaluator<FuncValueEvalTraits> FuncValueEval;

        struct FuncGradEvalTraits
        {
          enum
          {
            image_dim = TrafoEvalTraits::image_dim
          };
          typedef TrafoEvalType TrafoEvaluator;
          typedef TrafoEvalData TrafoData;
          typedef DataType_ DataType;
          typedef Tiny::Vector<DataType_, image_dim> ValueType;
        };

        typedef typename Functor_::template GradientEvaluator<FuncGradEvalTraits> FuncGradEval;

        TrafoEvalType _trafo_eval;
        FuncValueEval _func_value;
        FuncGradEval _func_grad;
        DomainPointType _dom_point;

      public:
        explicit NodeFunctional(const Space_& space, const Functor_& functor) :
          BaseClass(space, functor),
          _trafo_eval(space.get_trafo()),
          _func_value(functor),
          _func_grad(functor)
        {
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);
          _func_value.prepare(_trafo_eval);
          _func_grad.prepare(_trafo_eval);
        }

        void finish()
        {
          _func_grad.finish();
          _func_value.finish();
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
          // map point
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, _dom_point);
          DataType_ value(DataType_(0));

          // evaluate value or gradient
          if(assign_idx == 0)
            _func_value(value, trafo_data);
          else
            _func_grad(value, trafo_data);

          // return value
          return value;
        }
      };

    } // namespace BognerFoxSchmit
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_BOGNER_FOX_SCHMIT_NODE_FUNCTIONAL_HPP
