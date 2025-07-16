// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/assembly/function_integral_info.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Assembly job for the integration of an analytic function
     *
     * \tparam DataType_
     * The (scalar) datatype in which the assembly is to be performed.
     *
     * \tparam Function_
     * The analytic function class that is to be integrated.
     *
     * \tparam Trafo_
     * The trafo on which to integrate.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Function_, typename Trafo_, int max_der_>
    class AnalyticFunctionIntegralJob
    {
    public:
      static_assert(Function_::can_value, "function cannot compute values");
      static_assert(Function_::can_grad || (max_der_ < 1), "function gradients are required but not available");
      static_assert(Function_::can_hess || (max_der_ < 2), "function hessians are required but not available");

      /// the datatype to be used by the assembly
      typedef DataType_ DataType;

      /// declare our analytic eval traits
      typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

      /// our function integral type
      typedef typename AnalyticFunctionIntegral<DataType, Function_>::Type FunctionIntegralType;

    public:
      /**
       * \brief Analytic Function Integral Task implementation
       */
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
        /// trafo evaluator
        typedef typename Trafo_::template Evaluator<typename Trafo_::ShapeType, DataType>::Type TrafoEvaluator;
        /// trafo evaluator
        TrafoEvaluator trafo_eval;
        /// trafo eval data type
        typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType trafo_data;
        /// the cubature rule
        typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType cubature_rule;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integrals
        FunctionIntegralType loc_integral;
        /// the accumulated job integral
        FunctionIntegralType& job_integral;

      public:
        explicit Task(AnalyticFunctionIntegralJob& job) :
          trafo_eval(job._trafo),
          trafo_data(),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          func_eval(job._function),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          trafo_eval.prepare(cell);
        }

        void assemble()
        {
          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // call helper to integrate
            Intern::AnaFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, func_eval, trafo_data.img_point);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the trafo to use
      const Trafo_& _trafo;
      /// the function to be integrated
      const Function_& _function;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the analytic function that is to be integrated.
       *
       * \param[in] trafo
       * A \resident reference to the trafo representing the domain over which to integrate.
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit AnalyticFunctionIntegralJob(const Function_& function, const Trafo_& trafo, const String& cubature) :
        _trafo(trafo),
        _function(function),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      /**
       * \brief Returns the assembled function integral info
       */
      FunctionIntegralType& result()
      {
        return _integral;
      }
    }; // class AnalyticFunctionIntegralJob<...>

    /**
     * \brief Assembly job for the integration of a discrete finite element function
     *
     * \tparam Vector_
     * The type of the coefficient vector of the discrete function. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element space used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Vector_, typename Space_, int max_der_>
    class DiscreteFunctionIntegralJob
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

    public:
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr SpaceTags space_tags = SpaceTags::value |
          (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
          (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);

        /// our assembly traits
        typedef Assembly::AsmTraits1<DataType, Space_, TrafoTags::jac_det, space_tags> AsmTraits;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the finite element space to be used
        const Space_& space;
        /// the cubature factory used for integration
        const typename AsmTraits::TrafoType& trafo;
        /// the trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the cubature rule used for integration
        typename AsmTraits::CubatureRuleType cubature_rule;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename Vector_::GatherAxpy gather_axpy;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(DiscreteFunctionIntegralJob& job) :
          vector(job._vector),
          space(job._space),
          trafo(space.get_trafo()),
          trafo_eval(trafo),
          space_eval(space),
          dof_mapping(space),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          trafo_data(),
          space_data(),
          local_vector(),
          gather_axpy(vector),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          // prepare dof mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);
        }

        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // do the dirty work
            Intern::DiscFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, space_data, local_vector, num_loc_dofs);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the vector that is to be integrated
      const Vector_& _vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector
       * A \resident reference to the coefficient vector of the discrete function.
       *
       * \param[in] space
       * A \resident reference to the finite element space used for the discretization
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit DiscreteFunctionIntegralJob(const Vector_& vector, const Space_& space, const String& cubature) :
        _vector(vector),
        _space(space),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      FunctionIntegralType& result()
      {
        return _integral;
      }
    }; // class DiscreteFunctionIntegralJob<...>

    /**
     * \brief Assembly job for the integration of a analytic vs discrete error function
     *
     * This class is used to integrate the error function (u - u_h), where u is an analytic function
     * and u_h is a discrete finite element function.
     *
     * \tparam Function_
     * The type of the analytic function u that is to be integrated.
     *
     * \tparam Vector_
     * The type of the coefficient vector of the discrete function u_h. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element space used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Function_, typename Vector_, typename Space_, int max_der_>
    class ErrorFunctionIntegralJob
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      /// make sure that function and vector are compatible
      static_assert(Intern::ErrCompatHelper<Function_, Vector_>::valid, "function and vector are incompatible");

      /// our function integral type
      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

      /// declare our analytic eval traits
      typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

    public:
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
        static constexpr SpaceTags space_tags = SpaceTags::value |
          (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
          (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);

        /// our assembly traits
        typedef Assembly::AsmTraits1<DataType, Space_, trafo_config, space_tags> AsmTraits;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the finite element space to be used
        const Space_& space;
        /// the cubature factory used for integration
        const typename AsmTraits::TrafoType& trafo;
        /// the trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the cubature rule used for integration
        typename AsmTraits::CubatureRuleType cubature_rule;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename Vector_::GatherAxpy gather_axpy;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(ErrorFunctionIntegralJob& job) :
          vector(job._vector),
          space(job._space),
          trafo(space.get_trafo()),
          trafo_eval(trafo),
          space_eval(space),
          dof_mapping(space),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          trafo_data(),
          space_data(),
          local_vector(),
          gather_axpy(vector),
          func_eval(job._function),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          // prepare dof mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);
        }

        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // do the dirty work
            Intern::ErrFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, func_eval, trafo_data.img_point,
              space_data, local_vector, num_loc_dofs);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the function to be integrated
      const Function_& _function;
      /// the vector that is to be integrated
      const Vector_& _vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the analytic function u.
       *
       * \param[in] vector
       * A \resident reference to the coefficient vector of the discrete function u_h.
       *
       * \param[in] space
       * A \resident reference to the finite element space used for the discretization of u_h
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit ErrorFunctionIntegralJob(const Function_& function, const Vector_& vector, const Space_& space, const String& cubature) :
        _function(function),
        _vector(vector),
        _space(space),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      FunctionIntegralType& result()
      {
        return _integral;
      }
    }; // class ErrorFunctionIntegralJob<...>

    namespace Intern
    {
      template<typename DT_, typename IT_, int max_der_>
      struct OutVectorHelper
      {
        typedef LAFEM::DenseVectorBlocked<DT_, IT_, max_der_+1> VectorType;
      };

      template<typename DT_, typename IT_>
      struct OutVectorHelper<DT_, IT_, 0>
      {
        typedef LAFEM::DenseVector<DT_, IT_> VectorType;
      };
    }

    /**
     * \brief Assembly job for the elementwise integration of a analytic vs discrete error function
     *
     * This class is used to integrate the error function (u - u_h), where u is an analytic function
     * and u_h is a discrete finite element function for each cell seperatly and saving this into
     * a corresponding vector.
     *
     * \tparam Function_
     * The type of the analytic function u that is to be integrated.
     *
     * \tparam Vector_
     * The type of the coefficient vector of the discrete function u_h. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element space used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Maximilian Esser
     */
    template<typename Function_, typename Vector_, typename Space_, int max_der_>
    class CellErrorFunctionIntegralJob
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      /// make sure that function and vector are compatible
      static_assert(Intern::ErrCompatHelper<Function_, Vector_>::valid, "function and vector are incompatible");

      /// our function integral type
      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

      /// declare our analytic eval traits
      typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

      /// our output vectortype
      typedef typename Intern::OutVectorHelper<DataType, Index, max_der_>::VectorType OutVectorType;

      /// our cell fucntion integral type
      typedef FunctionCellIntegralInfo<FunctionIntegralType, OutVectorType> FunctionCellIntegralType;

    public:
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = true;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
        static constexpr SpaceTags space_tags = SpaceTags::value |
          (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
          (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);

        /// our assembly traits
        typedef Assembly::AsmTraits1<DataType, Space_, trafo_config, space_tags> AsmTraits;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the finite element space to be used
        const Space_& space;
        /// the output vector
        OutVectorType& out_vector;
        /// the cubature factory used for integration
        const typename AsmTraits::TrafoType& trafo;
        /// the trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the cubature rule used for integration
        typename AsmTraits::CubatureRuleType cubature_rule;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename Vector_::GatherAxpy gather_axpy;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;
        /// the local integral
        FunctionIntegralType cell_integral;

      public:
        explicit Task(CellErrorFunctionIntegralJob& job) :
          vector(job._vector),
          space(job._space),
          out_vector(job._out_vec),
          trafo(space.get_trafo()),
          trafo_eval(trafo),
          space_eval(space),
          dof_mapping(space),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          trafo_data(),
          space_data(),
          local_vector(),
          gather_axpy(vector),
          func_eval(job._function),
          loc_integral(),
          job_integral(job._integral),
          cell_integral()
        {
        }

        void prepare(Index cell)
        {
          // prepare dof mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);
        }

        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // format cell local integral value
          cell_integral.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // do the dirty work
            Intern::ErrFunIntJobHelper<max_der_>::work(cell_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, func_eval, trafo_data.img_point,
              space_data, local_vector, num_loc_dofs);
          }
          // add our cell value
          loc_integral.push(cell_integral);
        }

        void scatter()
        {
          // get cell index
          Index cell = dof_mapping.get_current_cell_index();
          if constexpr(max_der_ == 0)
          {
            out_vector(cell, cell_integral.norm_h0_sqr);
          }
          else
          {
            // construct tiny vector
            Tiny::Vector<DataType, max_der_+1> loc_vec;
            loc_vec[0] = cell_integral.norm_h0_sqr;
            loc_vec[1] = cell_integral.norm_h1_sqr;
            if constexpr(max_der_ > 1)
              loc_vec[2] = cell_integral.norm_h2_sqr;
            //scatter the tiny vector
            out_vector(cell, loc_vec);
          }
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the function to be integrated
      const Function_& _function;
      /// the vector that is to be integrated
      const Vector_& _vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the vector of the cell data
      OutVectorType _out_vec;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A \resident reference to the analytic function u.
       *
       * \param[in] vector
       * A \resident reference to the coefficient vector of the discrete function u_h.
       *
       * \param[in] space
       * A \resident reference to the finite element space used for the discretization of u_h
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit CellErrorFunctionIntegralJob(const Function_& function, const Vector_& vector, const Space_& space, const String& cubature) :
        _function(function),
        _vector(vector),
        _space(space),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
        _out_vec = OutVectorType(_space.get_mesh().get_num_elements());
        _out_vec.format();
      }

      /// \returns A cellerror object containing all error data
      FunctionCellIntegralType result()
      {
        FunctionCellIntegralType output(_integral, std::move(_out_vec));
        _out_vec = OutVectorType(_space.get_mesh().get_num_elements());
        _out_vec.format();
        return output;
      }
    }; // class ErrorFunctionIntegralJob<...>
  } // namespace Assembly
} // namespace FEAT
