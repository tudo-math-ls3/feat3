// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/node_functional_base.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Hermite3
    {
      template<
        typename Space_,
        typename Shape_,
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
      class NodeFunctional<Space_, Shape::Simplex<2>, 0, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 3;

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

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          // compute trafo data
          DomainPointType dom_point(DataType_(0));
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          const auto func_value = func_eval.value(trafo_data.img_point);
          const auto func_grad = func_eval.gradient(trafo_data.img_point);

          // set node functional values
          node_data[0] = func_value;
          node_data[1] = func_grad[0];
          node_data[2] = func_grad[1];
        }
      };

      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Simplex<2>, 2, DataType_> :
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
        typedef typename TrafoType::template Evaluator<Shape::Simplex<2>, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          // compute trafo data
          DomainPointType dom_point(DataType_(1) / DataType_(3));
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data.img_point);
        }
      };

      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<1>, 0, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 2;

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

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          // compute trafo data
          DomainPointType dom_point(DataType_(0));
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          const auto func_value = func_eval.value(trafo_data.img_point);
          const auto func_grad = func_eval.gradient(trafo_data.img_point);

          // set node functional values
          node_data[0] = func_value;
          node_data[1] = func_grad[0];
        }
      };

      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<2>, 0, DataType_> :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = 4;

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

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

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
          DomainPointType dom_point(DataType_(0));
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          const auto func_value = func_eval.value(trafo_data.img_point);
          const auto func_grad = func_eval.gradient(trafo_data.img_point);
          const auto func_hess = func_eval.hessian(trafo_data.img_point);

          // set node functional values
          node_data[0] = func_value;
          node_data[1] = func_grad[0];
          node_data[2] = func_grad[1];
          node_data[3] = func_hess[0][1];
        }
      };
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAT
