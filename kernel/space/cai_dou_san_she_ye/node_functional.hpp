// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/node_functional_base.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace CaiDouSanSheYe
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
       * \brief Facet Node Functional implementation for Cai-Douglas-Santos-Sheen-Ye
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

        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

        TrafoEvalType _trafo_eval;
        DomainPointType _barycentre;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
        {
          _barycentre.format();
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
      * \brief Bubble Node Functional implementation for Cai-Douglas-Santos-Sheen-Ye
      *
      * \author Peter Zajac
      */
      template<
        typename Space_,
        int shape_dim_,
        typename DataType_>
      class NodeFunctional<Space_, Shape::Hypercube<shape_dim_>, 0, DataType_> :
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

        //static constexpr int image_dim = TrafoEvalTraits::image_dim;

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

          // declare our evaluation traits
          typedef Analytic::EvalTraits<DataType_, Function_> FuncEvalTraits;

          // declare function evaluator
          typename Function_::template Evaluator<FuncEvalTraits> func_eval(function);
          typename FuncEvalTraits::ValueType value(DataType_(0));

          TrafoEvalData trafo_data;

          // 2-point Gauss quadrature point
          const DataType g2 = DataType(FEAT_F128C(0.57735026918962576450914878050195745564760175127));

          DomainPointType dom_point;

          // integrate f(x,y)*x*y on the reference element
          for(int i(0); i < (1 << shape_dim_); ++i)
          {
            // set cubature shape_dim_
            DataType_ aux = DataType(1);
            for(int j(0); j < shape_dim_; ++j)
              aux *= (dom_point[j] = Shape::ReferenceCell<ShapeType>::template vertex<DataType>(i, j) * g2);
            _trafo_eval(trafo_data, dom_point);
            value += aux * func_eval.value(trafo_data.img_point);
          }

          // set integral mean
          node_data[0] = (DataType_(1) / DataType_(1 << shape_dim_)) * value;
        }
      };
    } // namespace CaiDouSanSheYe
  } // namespace Space
} // namespace FEAT
