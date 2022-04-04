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
#include <kernel/util/tiny_algebra.hpp>


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

        // The coef Matrix gets calculatet in prepare.
       // Every row represents the coefficients of a basis function.
       // The i-th basis functions can be evaluated in (x,y) via Intern::base_q(x,y,coef[i])
        Tiny::Matrix<DataType_, 3, 3> coef;

        static constexpr TrafoTags trafo_tags = TrafoTags::img_point | TrafoTags::jac_det;
        typedef typename TrafoEvalType::template ConfigTraits<trafo_tags>::EvalDataType TrafoEvalData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"))
        {
          coef.format();
        }

        // prepare() is called for every FE.
        // In prepare the coefficientmatrix is calculated for every FE because otherwise the isoparametic trafo doesnt work.
        void prepare(Index cell_index)
        {
          coef.format();
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);

          // declare trafo evaluation data
          TrafoEvalData trafo_data;


          Tiny::Matrix<DataType_, 3, 3> temp(DataType_(0));
          DataType mean(DataType(0));

          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            // get current cubature point
            _trafo_eval(trafo_data, _cub_rule.get_point(i));

            // calculate p0,p1,p2,q0,q1,q2 at the i-th cubature point[x] and temporary store them in temp like this:
            //temp=[p0(x),q0(x),0;
            //      p1(x),q1(x),0;
            //      p2(x),q2(x),0;]
            for (int j = 0; j < 3; ++j)
            {
              temp(j, 0) = Intern::p(j, trafo_data.dom_point[0]);
              temp(j, 1) = Intern::q(j, trafo_data.dom_point[0]);
            }

            DataType_ mult = _cub_rule.get_weight(i) * trafo_data.jac_det;
            mean += mult;
            for (int px = 0; px < 3; ++px)
            {
              for (int qx = 0; qx < 3; ++qx)
              {
                coef(px, qx) +=
                  (mult * temp(px, 0) * temp(qx, 1));
              } // end for qx
            } // end for px
            temp.format();
          } // end cubature
          // scale coefficient matrix
          coef = coef * (2 / mean);
          temp.set_inverse(coef);
          coef.set_transpose(temp);
        } // end prepare

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

          // integrate over facet with cubature factory
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point)
              *Intern::base_q(trafo_data.dom_point[0], coef[2]);
            mean += _cub_rule.get_weight(i) *trafo_data.jac_det;
          }

          // set node_data to the integral mean
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

        // The coef Matrix gets calculatet in prepare.
        // Every row represents the coefficients of a basis function.
        // The i-th basis functions can be evaluated in (x,y) via Intern::base_q(x,y,coef[i])
        Tiny::Matrix<DataType_, 9, 9> coef;


        // declare trafo evaluation data
        static constexpr TrafoTags trafo_tags = TrafoTags::img_point | TrafoTags::jac_det;
        typedef typename TrafoEvalType::template ConfigTraits<trafo_tags>::EvalDataType TrafoEvalData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"))
        {
          coef.format();
        }

        // prepare() is called for every FE.
        // In prepare the coefficientmatrix is calculated for every FE because of the changing jacobi determinant.
        void prepare(Index cell_index)
        {
          coef.format();
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);

          // declare trafo evaluation data
          TrafoEvalData trafo_data;

          Tiny::Matrix<DataType_, 9, 9> temp(DataType_(0));
          DataType mean(DataType(0));
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            // get current cubature point
            _trafo_eval(trafo_data, _cub_rule.get_point(i));

            // calculate p0,p1,p2,q0,q1,q2 at the i-th cubature point[x,y] and temporary store them in temp like this:
            //temp=[p0(x),p0(y),q0(x),q0(y),0,..0;
            //      p1(x),p1(y),q1(x),q1(y),0,..0;
            //      p2(x),p2(y),q2(x),q2(y),0,..0;
            //      0,......,0]
            for (int j = 0; j < 3; ++j)
            {
              temp(j, 0) = Intern::p(j, trafo_data.dom_point[0]);
              temp(j, 1) = Intern::p(j, trafo_data.dom_point[1]);
              temp(j, 2) = Intern::q(j, trafo_data.dom_point[0]);
              temp(j, 3) = Intern::q(j, trafo_data.dom_point[1]);
            }

            DataType_ mult = _cub_rule.get_weight(i) * trafo_data.jac_det;
            mean += mult;
            for (int px = 0; px < 3; ++px)
            {
              for (int py = 0; py < 3; ++py)
              {
                for (int qx = 0; qx < 3; ++qx)
                {
                  for (int qy = 0; qy < 3; ++qy)
                  {
                    coef(px * 3 + py, qx * 3 + qy) +=
                      (mult * temp(px, 0) * temp(py, 1) * temp(qx, 2) * temp(qy, 3));
                  } // end for qy
                } // end for qx
              } // end for p_y
            } // end for px
            temp.format();
          } // end cubature
          // scale coefficient matrix
          coef = coef * (4 / mean);
          temp.set_inverse(coef);
          coef.set_transpose(temp);
        } // end prepare

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


          // integrate over facet
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point)
              * Intern::base_q(trafo_data.dom_point[0], trafo_data.dom_point[1], coef[8]);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set node_date to integral mean
          node_data[0] = (DataType_(4) / mean) * value;
        }
      };

      template<
        typename Space_,
        typename DataType_>
        class NodeFunctional<Space_, Shape::Hypercube<3>, DataType_> :
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
        typedef Shape::Hypercube<3> ShapeType;
        typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvalType;
        typedef typename TrafoEvalType::EvalTraits TrafoEvalTraits;
        typedef typename TrafoEvalTraits::DataType DataType;
        typedef typename TrafoEvalTraits::DomainPointType DomainPointType;

        static constexpr int image_dim = TrafoEvalTraits::image_dim;

        // The coef Matrix gets calculatet in prepare.
        // Every row represents the coefficients of a basis function.
        // The i-th basis functions can be evaluated in (x,y) via Intern::base_q(x,y,coef[i])
        Tiny::Matrix<DataType_, 27, 27> coef;


        // declare trafo evaluation data
        static constexpr TrafoTags trafo_tags = TrafoTags::img_point | TrafoTags::jac_det;
        typedef typename TrafoEvalType::template ConfigTraits<trafo_tags>::EvalDataType TrafoEvalData;

        typedef Cubature::Rule<ShapeType, DataType_, DataType_, DomainPointType> CubRuleType;

        TrafoEvalType _trafo_eval;
        CubRuleType _cub_rule;

      public:
        explicit NodeFunctional(const Space_& space) :
          BaseClass(space),
          _trafo_eval(space.get_trafo()),
          _cub_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"))
        {
          coef.format();
        }

        // prepare() is called for every FE.
        // In prepare the coefficientmatrix is calculated for every FE because of the changing jacobi determinant.
        void prepare(Index cell_index)
        {
          coef.format();
          BaseClass::prepare(cell_index);
          _trafo_eval.prepare(cell_index);

          // declare trafo evaluation data
          TrafoEvalData trafo_data;

          Tiny::Matrix<DataType_, 27, 27> temp(DataType_(0));
          DataType mean(DataType(0));
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            // get current cubature point
            _trafo_eval(trafo_data, _cub_rule.get_point(i));

            // calculate p0,p1,p2,q0,q1,q2 at the i-th cubature point[x,y,z] and temporary store them in temp like this:
            //temp=[p0(x),p0(y),p0(z),q0(x),q0(y),q0(z),0,..0;
            //      p1(x),p1(y),p1(z),q1(x),q1(y),q1(z),0,..0;
            //      p2(x),p2(y),p2(z),q2(x),q2(y),q2(z),0,..0;
            //      0,......,0]
            for (int j = 0; j < 3; ++j)
            {
              temp(j, 0) = Intern::p(j, trafo_data.dom_point[0]);
              temp(j, 1) = Intern::p(j, trafo_data.dom_point[1]);
              temp(j, 2) = Intern::p(j, trafo_data.dom_point[2]);
              temp(j, 3) = Intern::q(j, trafo_data.dom_point[0]);
              temp(j, 4) = Intern::q(j, trafo_data.dom_point[1]);
              temp(j, 5) = Intern::q(j, trafo_data.dom_point[2]);
            }

            DataType_ mult = _cub_rule.get_weight(i) * trafo_data.jac_det;
            mean += mult;
            for (int pz = 0; pz < 3; ++pz)
            {
              for (int px = 0; px < 3; ++px)
              {
                for (int py = 0; py < 3; ++py)
                {
                  for (int qz = 0; qz < 3; ++qz)
                  {
                    for (int qx = 0; qx < 3; ++qx)
                    {
                      for (int qy = 0; qy < 3; ++qy)
                      {
                        coef(9*pz+px * 3 + py, 9*qz+qx * 3 + qy) +=
                          (mult * temp(px, 0) * temp(py, 1)*temp(pz,2) * temp(qx, 3) * temp(qy, 4)*temp(qz,5));
                      } // end for qy
                    } // end for qx
                  } // end for qz
                } // end for py
              } // end for px
            } //end for pz
            temp.format();
          } // end cubature
          // scale coefficient matrix
          coef = coef * (8 / mean);
          temp.set_inverse(coef);
          coef.set_transpose(temp);
        } // end prepare

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


          // integrate over facet
          for (int i(0); i < _cub_rule.get_num_points(); ++i)
          {
            _trafo_eval(trafo_data, _cub_rule.get_point(i));
            value += (_cub_rule.get_weight(i) * trafo_data.jac_det) * func_eval.value(trafo_data.img_point)
              * Intern::base_q(trafo_data.dom_point[0], trafo_data.dom_point[1], trafo_data.dom_point[2], coef[26]);
            mean += _cub_rule.get_weight(i) * trafo_data.jac_det;
          }

          // set node_date to integral mean
          node_data[0] = (DataType_(8) / mean) * value;
        }
      };
    } // namespace Bernstein2
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_BERNSTEIN2_NODE_FUNCTIONAL_HPP
