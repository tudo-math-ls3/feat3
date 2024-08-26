// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_Q1TBNP_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_Q1TBNP_NODE_FUNCTIONAL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/node_functional_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Q1TBNP
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
       * \brief Node Functional implementation for Q1TBNP
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

        template<typename Function_>
        struct Value
        {
          typedef typename Analytic::EvalTraits<DataType_, Function_>::ValueType Type;
        };

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Space_::ShapeType ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;
        typedef typename TrafoType::template Evaluator<FacetType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;

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
          static_assert(std::is_base_of<Analytic::Function, Function_>::value, "invalid function object");

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          typename FuncEvalTraits::ValueType value(DataType_(0));
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;

          // integrate over facet
          for(int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set integral mean
          node_data[0] = (DataType_(1) / mean) * value;
        }
      };

      /**
       * \brief Bubble Node Functional implementation for 2D Q1TBNP
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<2>, 0, DataType_> :
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
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;
        typedef typename TrafoEvalTraits::ImagePointType ImagePointType;
        typedef typename TrafoEvalTraits::JacobianInverseType JacobianInverseType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_inv>::EvalDataType InvLinTrafoData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

        JacobianInverseType _inv_lin_mat;
        ImagePointType _inv_lin_vec;

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

          DomainPointType dom_point(DataType(0));
          InvLinTrafoData trafo_data;
          _trafo_eval(trafo_data, dom_point);
          _inv_lin_mat = trafo_data.jac_inv;
          _inv_lin_vec = trafo_data.img_point;
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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          typename FuncEvalTraits::ValueType value(DataType_(0));
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;
          DomainPointType dom_point;

          // integrate over facet
          for(int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            dom_point.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point) * dom_point[0] * dom_point[1];
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set integral mean
          node_data[0] = (DataType_(1) / mean) * value;
        }
      };

      /**
       * \brief Bubble Node Functional implementation for 3D Q1TBNP
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<3>, 0, DataType_> :
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
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;
        typedef typename TrafoEvalTraits::ImagePointType ImagePointType;
        typedef typename TrafoEvalTraits::JacobianInverseType JacobianInverseType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_inv>::EvalDataType InvLinTrafoData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

        JacobianInverseType _inv_lin_mat;
        ImagePointType _inv_lin_vec;

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

          DomainPointType dom_point(DataType(0));
          InvLinTrafoData trafo_data;
          _trafo_eval(trafo_data, dom_point);
          _inv_lin_mat = trafo_data.jac_inv;
          _inv_lin_vec = trafo_data.img_point;
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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          Tiny::Vector<typename FuncEvalTraits::ValueType, 3> vals;
          vals.format();
          DataType_ mean(DataType_(0));
          TrafoEvalData trafo_data;
          DomainPointType dom_point;

          // integrate over facet
          for(int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            dom_point.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);
            typename FuncEvalTraits::ValueType fv = (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point);
            vals[0] += fv * (dom_point[0] * dom_point[1]);
            vals[1] += fv * (dom_point[1] * dom_point[2]);
            vals[2] += fv * (dom_point[2] * dom_point[0]);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set integral mean
          node_data[0] = (DataType_(1) / mean) * vals[0];
          node_data[1] = (DataType_(1) / mean) * vals[1];
          node_data[2] = (DataType_(1) / mean) * vals[2];
        }
      };
    } // namespace Q1TBNP
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_Q1TBNP_NODE_FUNCTIONAL_HPP
