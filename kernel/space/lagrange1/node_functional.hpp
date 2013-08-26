#pragma once
#ifndef KERNEL_SPACE_LAGRANGE1_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_LAGRANGE1_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Lagrange1
    {
      template<
        typename Space_,
        typename Functor_,
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
      class NodeFunctional<Space_, Functor_, 0, DataType_> :
        public NodeFunctionalBase<Space_, Functor_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Functor_, DataType_> BaseClass;

      protected:
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_img_point = 1
          };
        };

        typedef typename Space_::TrafoType TrafoType;
        typedef typename TrafoType::template Evaluator<Shape::Vertex, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        typedef Trafo::EvalData<TrafoEvalTraits, TrafoConfig> TrafoEvalData;

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

        TrafoEvalType _trafo_eval;
        FuncEval _func_eval;
        DomainPointType _dom_point;

      public:
        explicit NodeFunctional(const Space_& space, const Functor_& functor) :
          BaseClass(space, functor),
          _trafo_eval(space.get_trafo()),
          _func_eval(functor)
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
          TrafoEvalData trafo_data;
          trafo_data(_trafo_eval, _dom_point);
          DataType_ value(DataType_(0));
          _func_eval(value, trafo_data);
          return value;
        }
      };
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_LAGRANGE1_NODE_FUNCTIONAL_HPP
