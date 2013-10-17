#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_DISCONTINUOUS_NODE_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/space/node_functional_base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Discontinuous
    {
      /// \cond internal
      namespace Intern
      {
        template<typename Shape_>
        struct Barycentre;

        template<int dim_>
        struct Barycentre< Shape::Simplex<dim_> >
        {
          template<typename T_>
          static void make(T_& p)
          {
            typedef typename T_::DataType DT;
            for(Index i(0); i < Index(dim_); ++i)
              p[i] = DT(1) / DT(dim_+1);
          }
        };

        template<int dim_>
        struct Barycentre< Shape::Hypercube<dim_> >
        {
          template<typename T_>
          static void make(T_& p)
          {
            typedef typename T_::DataType DT;
            for(Index i(0); i < Index(dim_); ++i)
              p[i] = DT(0);
          }
        };
      } // namespace Intern
      /// \endcond

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
      class NodeFunctional<Space_, Functor_, 0, Variant_, DataType_> :
        public NodeFunctionalBase<Space_, Functor_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, Functor_, DataType_> BaseClass;

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
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

        //typedef typename TrafoEvalType::template TrafoConfig<TrafoConfig1> TrafoConfig;
        //typedef Trafo::EvalData<TrafoEvalTraits, TrafoConfig> TrafoEvalData;
        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

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
          // set cell midpoint
          Intern::Barycentre<ShapeType>::make(_dom_point);
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
          _trafo_eval(trafo_data, _dom_point);
          DataType_ value(DataType_(0));
          _func_eval(value, trafo_data);
          return value;
        }
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_NODE_FUNCTIONAL_HPP
