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
            for(int i(0); i < dim_; ++i)
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
            for(int i(0); i < dim_; ++i)
              p[i] = DT(0);
          }
        };
      } // namespace Intern
      /// \endcond

      template<
        typename Space_,
        int codim_,
        typename Variant_,
        typename DataType_,
        typename Shape_ = typename Space_::ShapeType>
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
        typename DataType_,
        typename Shape_>
      class NodeFunctional<Space_, 0, Variant::StdPolyP<0>, DataType_, Shape_> :
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
        typedef typename TrafoType::template Evaluator<Shape_, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        // declare trafo evaluation data
        typedef typename TrafoEvalType::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;

        TrafoEvalType _trafo_eval;
        DomainPointType _dom_point;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo())
        {
          // set cell midpoint
          Intern::Barycentre<Shape_>::make(_dom_point);
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
          TrafoEvalData trafo_data;
          _trafo_eval(trafo_data, _dom_point);

          // evaluate function
          node_data[0] = func_eval.value(trafo_data.img_point);
        }
      };

      template<
        typename Space_,
        typename DataType_,
        int shape_dim_>
      class NodeFunctional<Space_, 0, Variant::StdPolyP<1>, DataType_, Shape::Simplex<shape_dim_> > :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = shape_dim_+1;

        template<typename Function_>
        struct Value
        {
          typedef typename Analytic::EvalTraits<DataType_, Function_>::ValueType Type;
        };

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef Shape::Simplex<shape_dim_> ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
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
          TrafoEvalData trafo_data;

          // loop over all nodal points
          DomainPointType dom_point;
          for(int i(0); i < max_assigned_dofs; ++i)
          {
            // set up domain point
            for(int j(0); j < shape_dim_; ++j)
            {
              dom_point[j] = DataType_(j+1 == i ? 1 : 0);
            }

            // evaluate trafo
            _trafo_eval(trafo_data, dom_point);

            // evaluate function
            node_data[i] = func_eval.value(trafo_data.img_point);
          }
        }
      };

      template<
        typename Space_,
        typename DataType_,
        int shape_dim_>
      class NodeFunctional<Space_, 0, Variant::StdPolyP<1>, DataType_, Shape::Hypercube<shape_dim_> > :
        public NodeFunctionalBase<Space_, DataType_>
      {
      public:
        typedef NodeFunctionalBase<Space_, DataType_> BaseClass;
        static constexpr int max_assigned_dofs = shape_dim_+1;

        template<typename Function_>
        struct Value
        {
          typedef typename Analytic::EvalTraits<DataType_, Function_>::ValueType Type;
        };

      protected:
        typedef typename Space_::TrafoType TrafoType;
        typedef Shape::Hypercube<shape_dim_> ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
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
          TrafoEvalData trafo_data;

          // evaluate center point
          DomainPointType dom_point(DataType_(0));
          _trafo_eval(trafo_data, dom_point);
          node_data[0] = func_eval.value(trafo_data.img_point);

          // loop over all nodal points
          for(int i(0); i < shape_dim_; ++i)
          {
            // evaluate domain point #1
            dom_point[i] = DataType_(1);
            _trafo_eval(trafo_data, dom_point);
            const auto value_0 = func_eval.value(trafo_data.img_point);

            // evaluate domain point #2
            dom_point[i] = -DataType_(1);
            _trafo_eval(trafo_data, dom_point);
            const auto value_1 = func_eval.value(trafo_data.img_point);

            // scale
            node_data[i+1] = DataType_(0.5) * (value_0 - value_1);

            // reset domain point
            dom_point[i] = DataType_(0);
          }
        }
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAT
