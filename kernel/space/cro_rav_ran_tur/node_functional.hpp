// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP
#define KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/node_functional_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

namespace FEAT
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

        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

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
          static_assert(std::is_base_of<Analytic::Function, Function_>::value, "invalid function object");

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);

          // compute trafo data
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, _barycentre);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data.img_point);
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
        /// \compilerhack PGI and Intel(18) do not understand complex template statements
#if defined(FEAT_COMPILER_PGI) || (defined(FEAT_COMPILER_INTEL) && FEAT_COMPILER_INTEL >= 1800 && FEAT_COMPILER_INTEL < 1900)
        static constexpr TrafoTags trafo_tags = TrafoTags::img_point | TrafoTags::jac_det;
        typedef typename TrafoEvalType::template ConfigTraits<trafo_tags>::EvalDataType TrafoEvalData;
#else
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point | TrafoTags::jac_det>::EvalDataType TrafoEvalData;
#endif

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
    } // namespace CroRavRanTur
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_CRO_RAV_RAN_TUR_NODE_FUNCTIONAL_HPP
