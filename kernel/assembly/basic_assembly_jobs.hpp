// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_BASIC_ASSEMBLY_JOBS_HPP
#define KERNEL_ASSEMBLY_BASIC_ASSEMBLY_JOBS_HPP 1

#include <kernel/assembly/base.hpp>
#include <kernel/analytic/function.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Basic Vector assembly task CRTP base-class
     *
     * This CRTP class can be used as a base-class for various types of vector assemblies tasks.
     * This class implements all the five interface member functions required by the DomainAssemblyJob
     * interface. The #assemble member function calls the two functions #set_point and #eval via
     * static down-cast, so these two functions have to be provided by the derived class.
     *
     * \tparam Derived_
     * The most derived CRTP class
     *
     * \tparam Vector_
     * The type of the vector that is to be assembled
     *
     * \tparam Space_
     * The finite element space to be used as test-space
     *
     * \tparam trafo_config_
     * The trafo configuration tags required for the assembly
     *
     * \tparam space_config_
     * The space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename Vector_, typename Space_,
      TrafoTags trafo_config_, SpaceTags space_config_>
    class BasicVectorAssemblyTaskCRTP
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

    protected:
      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_ | TrafoTags::jac_det,
        space_config_
      > AsmTraits;

      /// the vector that is to be assembled
      Vector_& vector;
      /// the test-/trial-space to be used
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
      /// the vector scatter object
      typename Vector_::ScatterAxpy scatter_axpy;
      /// the scatter scaling factor
      DataType scatter_alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[inout] vector_
       * A \resident reference to the vector that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the finite element space to be used for the assembly
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit BasicVectorAssemblyTaskCRTP(Vector_& vector_, const Space_& space_,
        const Cubature::DynamicFactory& cubature_factory, DataType alpha_) :
        vector(vector_),
        space(space_),
        trafo(space.get_trafo()),
        trafo_eval(trafo),
        space_eval(space),
        dof_mapping(space),
        cubature_rule(Cubature::ctor_factory, cubature_factory),
        scatter_axpy(vector),
        scatter_alpha(alpha_)
      {
      }

#ifdef DOXYGEN
      /**
       * \brief Sets the current cubature point
       *
       * \attention
       * This function must be implemented in the derived CRTP class!
       *
       * \param[in] tau
       * The \transient trafo evaluation data for the current cubature point.
       */
      void set_point(typename AsmTraits::TrafoEvalData& tau);

      /**
       * \brief Evaluates the linearform for a test-basis function
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] psi
       * The \transient evaluation data of the current test basis function.
       */
      void eval(ValueType& val, const DataType& weight, const typename AsmTraits::SpaceBasisData& psi);
#endif // DOXYGEN

      /**
       * \brief Prepares the task for the assembly on a given mesh cell.
       *
       * \param[in] cell
       * The index of the cell that is to be assembled on
       */
      void prepare(Index cell)
      {
        // prepare dof mapping
        dof_mapping.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);
      }

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local vector
        local_vector.format();

        // fetch number of local dofs
        const int num_loc_dofs = space_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          space_eval(space_data, trafo_data);

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(trafo_data);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // evaluate functional and integrate
            static_cast<Derived_&>(*this).eval(local_vector(i),
              trafo_data.jac_det * cubature_rule.get_weight(k),
              space_data.phi[i]);
            // continue with next test function
          }
          // continue with next cubature point
        }
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        scatter_axpy(local_vector, dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finishes the assembly on the current cell.
       */
      void finish()
      {
        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();

        // finish dof mapping
        dof_mapping.finish();
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class BasicVectorAssemblyTaskCRTP<...>

    /**
     * \brief Basic Matrix assembly task CRTP base-class for identical test-/trial-space
     *
     * This CRTP class can be used as a base-class for various types of matrix assemblies tasks.
     * This class implements all the five interface member functions required by the DomainAssemblyJob
     * interface. The #assemble member function calls the two functions #set_point and #eval via
     * static down-cast, so these two functions have to be provided by the derived class.
     *
     * \tparam Derived_
     * The most derived CRTP class
     *
     * \tparam Matrix_
     * The type of the matrix that is to be assembled
     *
     * \tparam Space_
     * The finite element space to be used as test-/trial-space
     *
     * \tparam trafo_config_
     * The trafo configuration tags required for the assembly
     *
     * \tparam space_config_
     * The space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename Matrix_, typename Space_,
      TrafoTags trafo_config_, SpaceTags space_config_>
    class BasicMatrixAssemblyTaskCRTP1
    {
    public:
      /// the data-type of the matrix
      typedef typename Matrix_::DataType DataType;
      /// the value-type of the matrix
      typedef typename Matrix_::ValueType ValueType;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

    protected:
      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_ | TrafoTags::jac_det,
        space_config_
      > AsmTraits;

      /// the matrix that is to be assembled
      Matrix_& matrix;
      /// the test-/trial-space to be used
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
      /// the local matrix to be assembled
      typename AsmTraits::template TLocalMatrix<ValueType> local_matrix;
      /// the matrix scatter object
      typename Matrix_::ScatterAxpy scatter_axpy;
      /// the scatter scaling factor
      DataType scatter_alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the finite element space to be used for the assembly
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit BasicMatrixAssemblyTaskCRTP1(Matrix_& matrix_, const Space_& space_,
        const Cubature::DynamicFactory& cubature_factory, DataType alpha_) :
        matrix(matrix_),
        space(space_),
        trafo(space.get_trafo()),
        trafo_eval(trafo),
        space_eval(space),
        dof_mapping(space),
        cubature_rule(Cubature::ctor_factory, cubature_factory),
        scatter_axpy(matrix),
        scatter_alpha(alpha_)
      {
      }

#ifdef DOXYGEN
      /**
       * \brief Sets the current cubature point
       *
       * \attention
       * This function must be implemented in the derived CRTP class!
       *
       * \param[in] tau
       * The \transient trafo evaluation data for the current cubature point.
       */
      void set_point(typename AsmTraits::TrafoEvalData& tau);

      /**
       * \brief Evaluates the bilinearform for a test-/trial-basis function pair
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] phi
       * The \transient evaluation data of the current trial basis function.
       *
       * \param[in] psi
       * The \transient evaluation data of the current test basis function.
       */
      void eval(ValueType& val, const DataType& weight,
        const typename AsmTraits::TrialBasisData& phi, const typename AsmTraits::TestBasisData& psi);
#endif // DOXYGEN

      /**
       * \brief Prepares the task for the assembly on a given mesh cell.
       *
       * \param[in] cell
       * The index of the cell that is to be assembled on
       */
      void prepare(Index cell)
      {
        // prepare dof mapping
        dof_mapping.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);
      }

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local matrix
        local_matrix.format();

        // fetch number of local dofs
        const int num_loc_dofs = space_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          space_eval(space_data, trafo_data);

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(trafo_data);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_dofs; ++j)
            {
              // evaluate operator and integrate
              static_cast<Derived_&>(*this).eval(local_matrix(i,j),
                trafo_data.jac_det * cubature_rule.get_weight(k),
                space_data.phi[j], space_data.phi[i]);
              // continue with next trial function
            }
            // continue with next test function
          }
          // continue with next cubature point
        }
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        // incorporate local matrix
        scatter_axpy(local_matrix, dof_mapping, dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finishes the assembly on the current cell.
       */
      void finish()
      {
        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();

        // finish dof mapping
        dof_mapping.finish();
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class BasicMatrixAssemblyTaskCRTP1<...>

    /**
     * \brief Basic Matrix assembly task CRTP base-class for (possibly different) test-/trial-spaces
     *
     * This CRTP class can be used as a base-class for various types of matrix assemblies tasks.
     * This class implements all the five interface member functions required by the DomainAssemblyJob
     * interface. The #assemble member function calls the two functions #set_point and #eval via
     * static down-cast, so these two functions have to be provided by the derived class.
     *
     * \tparam Derived_
     * The most derived CRTP class
     *
     * \tparam Matrix_
     * The type of the matrix that is to be assembled
     *
     * \tparam TestSpace_
     * The finite element space to be used as test-space
     *
     * \tparam TrialSpace_
     * The finite element space to be used as trial-space
     *
     * \tparam trafo_config_
     * The trafo configuration tags required for the assembly
     *
     * \tparam test_config_
     * The test-space configuration tags required for the assembly
     *
     * \tparam trial_config_
     * The trial-space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<
      typename Derived_,
      typename Matrix_,
      typename TestSpace_,
      typename TrialSpace_,
      TrafoTags trafo_config_,
      SpaceTags test_config_,
      SpaceTags trial_config_>
    class BasicMatrixAssemblyTaskCRTP2
    {
    public:
      /// the data-type of the matrix
      typedef typename Matrix_::DataType DataType;
      /// the value-type of the matrix
      typedef typename Matrix_::ValueType ValueType;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

    protected:
      /// our assembly traits
      typedef Assembly::AsmTraits2<
        DataType,
        TestSpace_,
        TrialSpace_,
        trafo_config_ | TrafoTags::jac_det,
        test_config_,
        trial_config_
      > AsmTraits;

      /// the matrix that is to be assembled
      Matrix_& matrix;
      /// the test-space to be used
      const TestSpace_& test_space;
      /// the trial-space to be used
      const TrialSpace_& trial_space;
      /// the cubature factory used for integration
      const typename AsmTraits::TrafoType& trafo;
      /// the trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval;
      /// the space evaluators
      typename AsmTraits::TestEvaluator test_eval;
      typename AsmTraits::TrialEvaluator trial_eval;
      /// the space dof-mappings
      typename AsmTraits::TestDofMapping test_dof_mapping;
      typename AsmTraits::TrialDofMapping trial_dof_mapping;
      /// the cubature rule used for integration
      typename AsmTraits::CubatureRuleType cubature_rule;
      /// the trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;
      /// the space evaluation data
      typename AsmTraits::TestEvalData test_data;
      typename AsmTraits::TrialEvalData trial_data;
      /// the local matrix to be assembled
      typename AsmTraits::template TLocalMatrix<ValueType> local_matrix;
      /// the matrix scatter object
      typename Matrix_::ScatterAxpy scatter_axpy;
      /// the scatter scaling factor
      DataType scatter_alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] test_space_
       * A \resident reference to the finite element test-space to be used for the assembly
       *
       * \param[in] trial_space_
       * A \resident reference to the finite element trial-space to be used for the assembly
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit BasicMatrixAssemblyTaskCRTP2(Matrix_& matrix_,
        const TestSpace_& test_space_, const TrialSpace_& trial_space_,
        const Cubature::DynamicFactory& cubature_factory, DataType alpha_) :
        matrix(matrix_),
        test_space(test_space_),
        trial_space(trial_space_),
        trafo(test_space.get_trafo()),
        trafo_eval(trafo),
        test_eval(test_space),
        trial_eval(trial_space),
        test_dof_mapping(test_space),
        trial_dof_mapping(trial_space),
        cubature_rule(Cubature::ctor_factory, cubature_factory),
        scatter_axpy(matrix),
        scatter_alpha(alpha_)
      {
      }

#ifdef DOXYGEN
      /**
       * \brief Sets the current cubature point
       *
       * \attention
       * This function must be implemented in the derived CRTP class!
       *
       * \param[in] tau
       * The \transient trafo evaluation data for the current cubature point.
       */
      void set_point(typename AsmTraits::TrafoEvalData& tau);

      /**
       * \brief Evaluates the bilinearform for a test-/trial-basis function pair
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] phi
       * The \transient evaluation data of the current trial basis function.
       *
       * \param[in] psi
       * The \transient evaluation data of the current test basis function.
       */
      void eval(ValueType& val, const DataType& weight,
        const typename AsmTraits::TrialBasisData& phi, const typename AsmTraits::TestBasisData& psi);
#endif // DOXYGEN

      /**
       * \brief Prepares the task for the assembly on a given mesh cell.
       *
       * \param[in] cell
       * The index of the cell that is to be assembled on
       */
      void prepare(Index cell)
      {
        // prepare dof mapping
        test_dof_mapping.prepare(cell);
        trial_dof_mapping.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluators
        test_eval.prepare(trafo_eval);
        trial_eval.prepare(trafo_eval);
      }

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local matrix
        local_matrix.format();

        // fetch number of local dofs
        const int num_loc_test_dofs = test_eval.get_num_local_dofs();
        const int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          test_eval(test_data, trafo_data);
          trial_eval(trial_data, trafo_data);

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(trafo_data);

          // test function loop
          for(int i(0); i < num_loc_test_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_trial_dofs; ++j)
            {
              // evaluate operator and integrate
              static_cast<Derived_&>(*this).eval(local_matrix(i,j),
                trafo_data.jac_det * cubature_rule.get_weight(k),
                trial_data.phi[j], test_data.phi[i]);
              // continue with next trial function
            }
            // continue with next test function
          }
          // continue with next cubature point
        }
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        // incorporate local matrix
        scatter_axpy(local_matrix, test_dof_mapping, trial_dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finishes the assembly on the current cell.
       */
      void finish()
      {
        // finish evaluators
        trial_eval.finish();
        test_eval.finish();
        trafo_eval.finish();

        // finish dof mapping
        trial_dof_mapping.finish();
        test_dof_mapping.finish();
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class BasicMatrixAssemblyTaskCRTP2<...>

    /**
     * \brief Vector assembly job for LinearFunctional implementations
     *
     * This class implements the DomainAssemblyJob interface to assemble a linearform that is
     * provided as a LinearFunctional implementation.
     *
     * \tparam LinearFunctional_
     * The linear functional that is to be assembled
     *
     * \tparam Vector_
     * The vector that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename LinearFunctional_, typename Vector_, typename Space_>
    class LinearFunctionalAssemblyJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = LinearFunctional_::trafo_config;
      static constexpr SpaceTags test_config = LinearFunctional_::test_config;

      class Task :
        public BasicVectorAssemblyTaskCRTP<Task, Vector_, Space_, trafo_config, test_config>
      {
      protected:
        /// our base-class typedef
        typedef BasicVectorAssemblyTaskCRTP<Task, Vector_, Space_, trafo_config, test_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename LinearFunctional_::template Evaluator<AsmTraits> func_eval;

      public:
        explicit Task(LinearFunctionalAssemblyJob& job) :
          BaseClass(job.vector, job.space, job.cubature_factory, job.alpha),
          func_eval(job.linear_functional)
        {
        }

        void prepare(Index cell)
        {
          BaseClass::prepare(cell);
          func_eval.prepare(this->trafo_eval);
        }

        void set_point(typename AsmTraits::TrafoEvalData& tau)
        {
          func_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight, typename AsmTraits::TestBasisData& psi)
        {
          Tiny::axpy(val, func_eval.eval(psi), weight);
        }

        void finish()
        {
          func_eval.finish();
          BaseClass::finish();
        }
      }; // class Task

    protected:
      /// a reference to the linear functional that is to be assembled
      const LinearFunctional_& linear_functional;
      /// a reference to the vector that is to be assembled
      Vector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the scaling factor for the assembly.
      DataType alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] linear_functional_
       * A \resident reference to the linear functional that is to be assembled.
       *
       * \param[inout] vector_
       * A \resident reference to the vector that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the (test) space to be used for the discretization.
       *
       * \param[in] cubature_
       * The name of the cubature rule that is to be used for integration.
       *
       * \param[in] alpha_
       * The scaling factor for the assembly.
       */
      explicit LinearFunctionalAssemblyJob(const LinearFunctional_& linear_functional_,
        Vector_& vector_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        linear_functional(linear_functional_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class LinearFunctionalAssemblyJob<...>

    /**
     * \brief Vector assembly job for a force functional represented by an analytic function
     *
     * This class implements an assembly job for a right-hand-side vector that is assembled from
     * a force functional. Effectively, this class performs the same assembly operation as the
     * LinearFunctionalAssemblyJob class in combination with the ForceFunctional linearform,
     * however, the caller does not have to create a ForceFunctional object but instead simply
     * passes the analytic function representing the right-hand-side function to this class.
     *
     * \tparam Function_
     * A template class implementing the AnalyticFunction interface representing the function \e f.
     *
     * \tparam Vector_
     * The vector that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename Function_, typename Vector_, typename Space_>
    class ForceFunctionalAssemblyJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
      static constexpr SpaceTags test_config = SpaceTags::value;

      class Task :
        public BasicVectorAssemblyTaskCRTP<Task, Vector_, Space_, trafo_config, test_config>
      {
      protected:
        /// our base-class typedef
        typedef BasicVectorAssemblyTaskCRTP<Task, Vector_, Space_, trafo_config, test_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// declare our analytic eval traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the function value in the current point
        typename AnalyticEvalTraits::ValueType func_value;

      public:
        explicit Task(ForceFunctionalAssemblyJob& job) :
          BaseClass(job.vector, job.space, job.cubature_factory, job.alpha),
          func_eval(job.function)
        {
        }

        void set_point(typename AsmTraits::TrafoEvalData& tau)
        {
          func_value = func_eval.value(tau.img_point);
        }

        void eval(ValueType& val, const DataType& weight, typename AsmTraits::TestBasisData& psi)
        {
          Tiny::axpy(val, func_value, weight * psi.value);
        }
      }; // class Task

    protected:
      /// a reference to the linear functional that is to be assembled
      const Function_& function;
      /// a reference to the vector that is to be assembled
      Vector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the scaling factor for the assembly.
      DataType alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function_
       * A \resident reference to the right-hand-side function that is to be assembled as a functional.
       *
       * \param[inout] vector_
       * A \resident reference to the vector that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the (test) space to be used for the discretization.
       *
       * \param[in] cubature_
       * The name of the cubature rule that is to be used for integration.
       *
       * \param[in] alpha_
       * The scaling factor for the assembly.
       */
      explicit ForceFunctionalAssemblyJob(const Function_& function_,
        Vector_& vector_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        function(function_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class ForceFunctionalAssemblyJob<...>

    /**
     * \brief Matrix assembly job for BilinearOperator and identical test-/trial-spaces
     *
     * This class implements the DomainAssemblyJob interface to assemble a bilinearform that is
     * provided as a BilinearOperator implementation.
     *
     * \tparam BilinearOperator_
     * The bilinear operator that is to be assembled
     *
     * \tparam Matrix_
     * The matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename BilinearOperator_, typename Matrix_, typename Space_>
    class BilinearOperatorMatrixAssemblyJob1
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags space_config = BilinearOperator_::test_config | BilinearOperator_::trial_config;

      class Task :
        public BasicMatrixAssemblyTaskCRTP1<Task, Matrix_, Space_, trafo_config, space_config>
      {
      protected:
        /// our base-class typedef
        typedef BasicMatrixAssemblyTaskCRTP1<Task, Matrix_, Space_, trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(BilinearOperatorMatrixAssemblyJob1& job) :
          BaseClass(job.matrix, job.space, job.cubature_factory, job.alpha),
          oper_eval(job.bilinear_operator)
        {
        }

        void prepare(Index cell)
        {
          BaseClass::prepare(cell);
          oper_eval.prepare(this->trafo_eval);
        }

        void set_point(typename AsmTraits::TrafoEvalData& tau)
        {
          oper_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight,
          typename AsmTraits::TrialBasisData& phi, typename AsmTraits::TestBasisData& psi)
        {
          val += weight * oper_eval.eval(phi, psi);
        }

        void finish()
        {
          oper_eval.finish();
          BaseClass::finish();
        }
      }; // class Task

    protected:
      /// a reference to the bilinear operator that is to be assembled
      const BilinearOperator_& bilinear_operator;
      /// a reference to the matrix that is to be assembled
      Matrix_& matrix;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the scaling factor for the assembly.
      DataType alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] bilinear_operator_
       * A \resident reference to the bilinear operator that is to be assembled.
       *
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the test/trial space to be used for the discretization.
       *
       * \param[in] cubature_
       * The name of the cubature rule that is to be used for integration.
       *
       * \param[in] alpha_
       * The scaling factor for the assembly.
       */
      explicit BilinearOperatorMatrixAssemblyJob1(const BilinearOperator_& bilinear_operator_,
        Matrix_& matrix_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        matrix(matrix_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class BilinearOperatorMatrixAssemblyJob1<...>

    /**
     * \brief Matrix assembly job for BilinearOperator and (possibly different) test-/trial-spaces
     *
     * This class implements the DomainAssemblyJob interface to assemble a bilinearform that is
     * provided as a BilinearOperator implementation.
     *
     * \tparam BilinearOperator_
     * The bilinear operator that is to be assembled
     *
     * \tparam Matrix_
     * The matrix that is to be assembled.
     *
     * \tparam TestSpace_
     * The finite element space to be used as test-space.
     *
     * \tparam TrialSpace_
     * The finite element space to be used as trial-space.
     *
     * \author Peter Zajac
     */
    template<typename BilinearOperator_, typename Matrix_, typename TestSpace_, typename TrialSpace_>
    class BilinearOperatorMatrixAssemblyJob2
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags test_config  = BilinearOperator_::test_config;
      static constexpr SpaceTags trial_config = BilinearOperator_::trial_config;

      class Task :
        public BasicMatrixAssemblyTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, test_config, trial_config>
      {
      protected:
        /// our base-class typedef
        typedef BasicMatrixAssemblyTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, test_config, trial_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(BilinearOperatorMatrixAssemblyJob2& job) :
          BaseClass(job.matrix, job.test_space, job.trial_space, job.cubature_factory, job.alpha),
          oper_eval(job.bilinear_operator)
        {
        }

        void prepare(Index cell)
        {
          BaseClass::prepare(cell);
          oper_eval.prepare(this->trafo_eval);
        }

        void set_point(typename AsmTraits::TrafoEvalData& tau)
        {
          oper_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight,
          typename AsmTraits::TrialBasisData& phi, typename AsmTraits::TestBasisData& psi)
        {
          Tiny::axpy(val, oper_eval.eval(phi, psi), weight);
        }

        void finish()
        {
          oper_eval.finish();
          BaseClass::finish();
        }
      }; // class Task

    protected:
      /// a reference to the bilinear operator that is to be assembled
      const BilinearOperator_& bilinear_operator;
      /// a reference to the matrix that is to be assembled
      Matrix_& matrix;
      /// a reference to the finite element space to be used as test-space
      const TestSpace_& test_space;
      /// a reference to the finite element space to be used as trial-space
      const TrialSpace_& trial_space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the scaling factor for the assembly.
      DataType alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] bilinear_operator_
       * A \resident reference to the bilinear operator that is to be assembled.
       *
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] test_space_
       * A \resident reference to the test space to be used for the discretization.
       *
       * \param[in] trial_space_
       * A \resident reference to the trial space to be used for the discretization.
       *
       * \param[in] cubature_
       * The name of the cubature rule that is to be used for integration.
       *
       * \param[in] alpha_
       * The scaling factor for the assembly.
       */
      explicit BilinearOperatorMatrixAssemblyJob2(const BilinearOperator_& bilinear_operator_,
        Matrix_& matrix_, const TestSpace_& test_space_, const TrialSpace_& trial_space_,
        String cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        matrix(matrix_),
        test_space(test_space_),
        trial_space(trial_space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class BilinearOperatorMatrixAssemblyJob2<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BASIC_ASSEMBLY_JOBS_HPP
