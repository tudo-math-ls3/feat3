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
        int dim_,
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
        typename DataType_>
      class NodeFunctional<Space_, 0, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 6;

        template<typename Function_>
        struct Value
        {
          typedef typename Analytic::EvalTraits<DataType_, Function_>::ValueType Type;
        };

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<Shape::Vertex, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          static constexpr bool need_img_point = true;
        };

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        TrafoEvalType _trafo_eval;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
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
          static_assert(std::is_base_of<Analytic::Function, Function_>::value, "invalid function object");
          static_assert(Function_::can_value, "function cannot compute values");
          static_assert(Function_::can_grad, "function cannot compute gradients");
          static_assert(Function_::can_hess, "function cannot compute hessians");

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          // compute trafo data
          DomainPointType dom_point;
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          typename FuncEvalTraits::ValueType func_value;
          typename FuncEvalTraits::GradientType func_grad;
          typename FuncEvalTraits::HessianType func_hess;
          func_eval.value(func_value, trafo_data.img_point);
          func_eval.gradient(func_grad, trafo_data.img_point);
          func_eval.hessian(func_hess, trafo_data.img_point);

          // set node functional values
          node_data[0] = func_value;
          node_data[1] = func_grad[0];
          node_data[2] = func_grad[1];
          node_data[3] = func_hess[0][0];
          node_data[4] = func_hess[1][1];
          node_data[5] = func_hess[0][1];
        }
      };

      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, 1, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 1;

        template<typename Function_>
        struct Value
        {
          typedef typename Analytic::EvalTraits<DataType_, Function_>::ValueType Type;
        };

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename TrafoType::template Evaluator<Shape::Simplex<1>, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          static constexpr bool need_img_point = true;
          static constexpr bool need_jac_mat = true;
        };

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoConfig>::EvalDataType TrafoEvalData;

        TrafoEvalType _trafo_eval;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
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
          static_assert(std::is_base_of<Analytic::Function, Function_>::value, "invalid function object");
          static_assert(Function_::can_grad, "function cannot compute gradients");

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          DomainPointType dom_point(DataType_(0.5));
          TrafoEvalData trafo_data;

          _trafo_eval(trafo_data, dom_point);

          // compute edge normal
          DataType_ dnx = +trafo_data.jac_mat(1,0);
          DataType_ dny = -trafo_data.jac_mat(0,0);
          DataType_ dnl = DataType_(1) / Math::sqrt(dnx*dnx + dny*dny);
          dnx *= dnl;
          dny *= dnl;

          // evaluate function gradient
          typename FuncEvalTraits::GradientType func_grad;
          func_eval.gradient(func_grad, trafo_data.img_point);

          // set edge normal derivative
          node_data[0] = dnx * func_grad[0] + dny * func_grad[1];
        }
      };
    } // namespace Argyris
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_ARGYRIS_NODE_FUNCTIONAL_HPP
