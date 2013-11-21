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
        int codim_,
        typename Variant_,
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

      template<
        typename Space_,
        typename Variant_,
        typename DataType_>
      class NodeFunctional<Space_, 0, Variant_, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr Index max_assigned_dofs = Index(1);

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
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
        DomainPointType _dom_point;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
        {
          // set cell midpoint
          Intern::Barycentre<ShapeType>::make(_dom_point);
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

        Index get_num_assigned_dofs() const
        {
          return max_assigned_dofs;
        }

        template<
          typename NodeData_,
          typename Function_>
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
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, _dom_point);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data);

          // finish function evaluator
          func_eval.finish();
        }
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_NODE_FUNCTIONAL_HPP
