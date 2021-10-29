// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_BERNSTEIN2_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_BERNSTEIN2_NODE_FUNCTIONAL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/node_functional_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/space/bernstein2/scalar_basis.hpp>


namespace FEAT
{
  namespace Space
  {
    namespace Bernstein2
    {
      template<
        typename Space_,
        typename ShapeType_,
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
      class NodeFunctional<Space_, Shape::Vertex, DataType_> :
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
        typedef typename TrafoType::template Evaluator<Shape::Vertex, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
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
          DomainPointType dom_point(DataType(0));
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, dom_point);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data.img_point);
        }
      };

      template<
        typename Space_,
        typename DataType_>
        class NodeFunctional<Space_, Shape::Hypercube<1>, DataType_> :
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
        typedef Shape::Hypercube<1> ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"))
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
          typename FuncEvalTraits::ValueType value(DataType_(0));
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;
          //DomainPointType dom_point(DataType_(0));
          //_trafo_eval(trafo_data, dom_point);


          // integrate over facet
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point)
              * Intern::q2(trafo_data.dom_point[0]);
              //*(DataType_(-7.5)*trafo_data.dom_point[0] * trafo_data.dom_point[0] +DataType_(3));
            mean += _cub_rule.get_weight(i) *trafo_data.jac_det;
          }

          //std::cout << "MeanTest:" << mean << std::endl;
          // set integral mean
          node_data[0] = (DataType_(2) / mean)* value;
        }
      };

      template<
        typename Space_,
        typename DataType_>
        class NodeFunctional<Space_, Shape::Hypercube<2>, DataType_> :
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
        typedef Shape::Hypercube<2> ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"))
        {
        }

        void prepare(Index cell_index)
        {
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);
          typename TrafoEvalType::template ConfigTraits<TrafoTags::jac_det>::EvalDataType trafo_data;
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
          typename FuncEvalTraits::ValueType value(DataType_(0));
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;
          //DomainPointType dom_point(DataType_(0));
          //_trafo_eval(trafo_data, dom_point);


          // integrate over facet
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point)
              * Intern::q2(trafo_data.dom_point[0]) * Intern::q2(trafo_data.dom_point[1]);
              //* (DataType_(-7.5) * trafo_data.dom_point[0] * trafo_data.dom_point[0] + DataType_(3))
              //* (DataType_(-7.5) * trafo_data.dom_point[1] * trafo_data.dom_point[1] + DataType_(3));
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

         // std::cout << "Mean:" << mean << std::endl;
          //std::cout << trafo_data.dom_point << std::endl;
          // set integral mean
          node_data[0] = (DataType_(4) / mean) * value;
        }
      };
    } // namespace Bernstein2
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_BERNSTEIN2_NODE_FUNCTIONAL_HPP
