// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/function_integral_info.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

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
    class DomainAssemblyBasicVectorTaskCRTP
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
      explicit DomainAssemblyBasicVectorTaskCRTP(Vector_& vector_, const Space_& space_,
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
    }; // class DomainAssemblyBasicVectorTaskCRTP<...>

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
    class DomainAssemblyBasicMatrixTaskCRTP1
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
      explicit DomainAssemblyBasicMatrixTaskCRTP1(Matrix_& matrix_, const Space_& space_,
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
    }; // class DomainAssemblyBasicMatrixTaskCRTP1<...>

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
    class DomainAssemblyBasicMatrixTaskCRTP2
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
      explicit DomainAssemblyBasicMatrixTaskCRTP2(Matrix_& matrix_,
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
    }; // class DomainAssemblyBasicMatrixTaskCRTP2<...>

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
    class DomainAssemblyLinearFunctionalVectorJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = LinearFunctional_::trafo_config;
      static constexpr SpaceTags test_config = LinearFunctional_::test_config;

      class Task :
        public DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, test_config>
      {
      protected:
        /// our base-class typedef
        typedef DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, test_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename LinearFunctional_::template Evaluator<AsmTraits> func_eval;

      public:
        explicit Task(DomainAssemblyLinearFunctionalVectorJob& job) :
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
      explicit DomainAssemblyLinearFunctionalVectorJob(const LinearFunctional_& linear_functional_,
        Vector_& vector_, const Space_& space_, const String& cubature_, DataType alpha_ = DataType(1)) :
        linear_functional(linear_functional_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class DomainAssemblyLinearFunctionalVectorJob<...>

    /**
     * \brief Assembles a linear functional into a vector
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[inout] vector
     * The \transient vector that is to be assembled.
     *
     * \param[in] linear_functional
     * A \transient reference to the functional implementing the Assembly::LinearFunctional interface to be assembled.
     *
     * \param[in] space
     * A \transient reference to the finite-element to be used as the test-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Vector_, typename LinFunc_, typename Space_>
    void assemble_linear_functional_vector(DomainAssembler<Trafo_>& dom_asm,
      Vector_& vector, const LinFunc_& linear_functional, const Space_& space,
      const String& cubature, const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector length");

      DomainAssemblyLinearFunctionalVectorJob<LinFunc_, Vector_, Space_> job(
        linear_functional, vector, space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Vector assembly job for a force functional represented by an analytic function
     *
     * This class implements an assembly job for a right-hand-side vector that is assembled from
     * a force functional. Effectively, this class performs the same assembly operation as the
     * DomainAssemblyLinearFunctionalVectorJob class in combination with the ForceFunctional linearform,
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
    class DomainAssemblyForceFunctionalVectorJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
      static constexpr SpaceTags test_config = SpaceTags::value;

      class Task :
        public DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, test_config>
      {
      protected:
        /// our base-class typedef
        typedef DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, test_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// declare our analytic eval traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the function value in the current point
        typename AnalyticEvalTraits::ValueType func_value;

      public:
        explicit Task(DomainAssemblyForceFunctionalVectorJob& job) :
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
      explicit DomainAssemblyForceFunctionalVectorJob(const Function_& function_,
        Vector_& vector_, const Space_& space_, const String& cubature_, DataType alpha_ = DataType(1)) :
        function(function_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class DomainAssemblyForceFunctionalVectorJob<...>

    /**
     * \brief Assembles a force function into a vector
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[inout] vector
     * The \transient vector that is to be assembled.
     *
     * \param[in] function
     * A \transient reference to the force function implementing the Analytic::Function interface to be assembled.
     *
     * \param[in] space
     * A \transient reference to the finite-element to be used as the test-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Vector_, typename Function_, typename Space_>
    void assemble_force_function_vector(DomainAssembler<Trafo_>& dom_asm,
      Vector_& vector, const Function_& function, const Space_& space,
      const String& cubature, const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector length");

      DomainAssemblyForceFunctionalVectorJob<Function_, Vector_, Space_> job(
        function, vector, space, cubature, alpha);
      dom_asm.assemble(job);
    }

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
    class DomainAssemblyBilinearOperatorMatrixJob1
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags space_config = BilinearOperator_::test_config | BilinearOperator_::trial_config;

      class Task :
        public DomainAssemblyBasicMatrixTaskCRTP1<Task, Matrix_, Space_, trafo_config, space_config>
      {
      protected:
        /// our base-class typedef
        typedef DomainAssemblyBasicMatrixTaskCRTP1<Task, Matrix_, Space_, trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(DomainAssemblyBilinearOperatorMatrixJob1& job) :
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
      explicit DomainAssemblyBilinearOperatorMatrixJob1(const BilinearOperator_& bilinear_operator_,
        Matrix_& matrix_, const Space_& space_, const String& cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        matrix(matrix_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class DomainAssemblyBilinearOperatorMatrixJob1<...>

    /**
     * \brief Assembles a bilinear operator into a matrix with identical test- and trial-spaces.
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[inout] matrix
     * The \transient matrix that is to be assembled.
     *
     * \param[in] bilinear_operator
     * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
     *
     * \param[in] space
     * A \transient reference to the finite-element to be used as the test- and trial-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Matrix_, typename BilOp_, typename Space_>
    void assemble_bilinear_operator_matrix_1(DomainAssembler<Trafo_>& dom_asm, Matrix_& matrix,
      const BilOp_& bilinear_operator, const Space_& space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix row count");

      DomainAssemblyBilinearOperatorMatrixJob1<BilOp_, Matrix_, Space_> job(
        bilinear_operator, matrix, space, cubature, alpha);
      dom_asm.assemble(job);
    }

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
    class DomainAssemblyBilinearOperatorMatrixJob2
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags test_config  = BilinearOperator_::test_config;
      static constexpr SpaceTags trial_config = BilinearOperator_::trial_config;

      class Task :
        public DomainAssemblyBasicMatrixTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, test_config, trial_config>
      {
      protected:
        /// our base-class typedef
        typedef DomainAssemblyBasicMatrixTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, test_config, trial_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(DomainAssemblyBilinearOperatorMatrixJob2& job) :
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
      explicit DomainAssemblyBilinearOperatorMatrixJob2(const BilinearOperator_& bilinear_operator_,
        Matrix_& matrix_, const TestSpace_& test_space_, const TrialSpace_& trial_space_,
        const String& cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        matrix(matrix_),
        test_space(test_space_),
        trial_space(trial_space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class DomainAssemblyBilinearOperatorMatrixJob2<...>

    /**
     * \brief Assembles a bilinear operator into a matrix with different test- and trial-spaces.
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[inout] matrix
     * The \transient matrix that is to be assembled.
     *
     * \param[in] bilinear_operator
     * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
     *
     * \param[in] test_space
     * A \transient reference to the finite-element test-space to be used.
     *
     * \param[in] trial_space
     * A \transient reference to the finite-element trial-space to be used.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Matrix_, typename BilOp_, typename TestSpace_, typename TrialSpace_>
    void assemble_bilinear_operator_matrix_2(DomainAssembler<Trafo_>& dom_asm,
      Matrix_& matrix, const BilOp_& bilinear_operator, const TestSpace_& test_space,
      const TrialSpace_& trial_space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == test_space.get_trafo(), "domain assembler and test space have different trafos");
      XASSERTM(dom_asm.get_trafo() == trial_space.get_trafo(), "domain assembler and trial space have different trafos");
      XASSERTM(matrix.columns() == trial_space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == test_space.get_num_dofs(), "invalid matrix row count");

      DomainAssemblyBilinearOperatorMatrixJob2<BilOp_, Matrix_, TestSpace_, TrialSpace_> job(
        bilinear_operator, matrix, test_space, trial_space, cubature, alpha);
      dom_asm.assemble(job);
    }

    /**
     * \brief Vector assembly job for BilinearOperator and identical test-/trial-spaces
     *
     * This class implements the DomainAssemblyJob interface to assemble a vector that is
     * calculated by
     *
     * \f$ \int_\Omega \tilde{v} \psi = \int_\Omega A(v, \psi) dx, m=1,\dots,d \f$
     * where \f$ A(\cdot, \cdot) \f$ is provided by a BilinearOperator interface and
     * \f$v\f$ is a finite element vector.
     *
     * \tparam BilinearOperator_
     * The bilinear operator that is to be assembled
     *
     * \tparam Vector_
     * The vector that is to be assembled.
     *
     * \tparam VectorSol_
     * The vector that the biliniear operator is to be applied to.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Maximilian Esser
     */
    template<typename BilinearOperator_, typename Vector_, typename VectorSol_, typename Space_>
    class DomainAssemblyBilinearOperatorApplyVectorJob1
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags space_config = BilinearOperator_::test_config | BilinearOperator_::trial_config;

      class Task :
        public DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, space_config>
      {
      protected:
        /// our base-class typedef
        typedef DomainAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;
        typedef typename BilinearOperator_::template Evaluator<AsmTraits>::ValueType MatValType;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

        /// the local solution vector that the bilinear operator is applied to
        const VectorSol_& vec_sol;

        /// the gather axpy for the local solution vector
        typename VectorSol_::GatherAxpy sol_gather;

        /// the local matrix that is used to assemble the local vector
        typename AsmTraits::template TLocalMatrix<MatValType> local_matrix;
        typename AsmTraits::template TLocalVector<ValueType> local_vec_sol;


      public:
        explicit Task(DomainAssemblyBilinearOperatorApplyVectorJob1& job) :
          BaseClass(job.vector, job.space, job.cubature_factory, job.alpha),
          oper_eval(job.bilinear_operator),
          vec_sol(job.vec_sol),
          sol_gather(vec_sol)
        {
        }

        void prepare(Index cell)
        {
          BaseClass::prepare(cell);
          oper_eval.prepare(this->trafo_eval);
          local_vec_sol.format();
          sol_gather(local_vec_sol, this->dof_mapping);
        }

        void set_point(typename AsmTraits::TrafoEvalData& tau)
        {
          oper_eval.set_point(tau);
        }

        void eval(MatValType& val, const DataType& weight,
          typename AsmTraits::TrialBasisData& phi, typename AsmTraits::TestBasisData& psi)
        {
          val += weight * oper_eval.eval(phi, psi);
        }

        /**
        * \brief Performs the local assembly.
        */
        void assemble()
        {
          // format local matrix
          local_matrix.format();
          this->local_vector.format();

          // fetch number of local dofs
          const int num_loc_dofs = this->space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < this->cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            this->trafo_eval(this->trafo_data, this->cubature_rule.get_point(k));

            // compute basis function data
            this->space_eval(this->space_data, this->trafo_data);

            // evaluate operator
            this->set_point(this->trafo_data);

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // evaluate operator and integrate
                this->eval(local_matrix(i,j),
                  this->trafo_data.jac_det * this->cubature_rule.get_weight(k),
                  this->space_data.phi[j], this->space_data.phi[i]);
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }
          // apply to local vector
          for(int i(0); i < num_loc_dofs; ++i)
            for(int j(0); j < num_loc_dofs; ++j)
              this->local_vector[i].add_mat_vec_mult(this->local_matrix(i,j), this->local_vec_sol[j]);
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
      /// a reference to the vector that is to be assembled
      Vector_& vector;
      /// a reference to the vector the biliniear operator is applied on
      const VectorSol_& vec_sol;
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
      explicit DomainAssemblyBilinearOperatorApplyVectorJob1(const BilinearOperator_& bilinear_operator_,
        Vector_& vector_, const VectorSol_& vec_sol_, const Space_& space_, const String& cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        vector(vector_),
        vec_sol(vec_sol_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class DomainAssemblyBilinearOperatorApplyVectorJob1<...>

    /**
     * \brief Assembles the application of a bilinear operator into a vector with identical test- and trial-spaces.
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[inout] vector
     * The \transient vector that is to be assembled.
     *
     * \param[in] vector_sol
     * The \transient vector that is to be applied upon.
     *
     * \param[in] bilinear_operator
     * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
     *
     * \param[in] space
     * A \transient reference to the finite-element to be used as the test- and trial-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] alpha
     * The scaling factor for the assembly.
     */
    template<typename Trafo_, typename Vector_, typename VectorSol_, typename BilOp_, typename Space_>
    void assemble_bilinear_operator_apply_vector_1(DomainAssembler<Trafo_>& dom_asm, Vector_& vector,
      const VectorSol_& vec_sol, const BilOp_& bilinear_operator, const Space_& space, const String& cubature,
      const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

      DomainAssemblyBilinearOperatorApplyVectorJob1<BilOp_, Vector_, VectorSol_, Space_> job(
        bilinear_operator, vector, vec_sol, space, cubature, alpha);
      dom_asm.assemble(job);
    }

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
    class DomainAssemblyAnalyticFunctionIntegralJob
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
        explicit Task(DomainAssemblyAnalyticFunctionIntegralJob& job) :
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
      explicit DomainAssemblyAnalyticFunctionIntegralJob(const Function_& function, const Trafo_& trafo, const String& cubature) :
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
    }; // class DomainAssemblyAnalyticFunctionIntegralJob<...>

    /**
     * \brief Assembles the integral of an analytic function
     *
     * \tparam DataType_
     * The (scalar) datatype in which the assembly is to be performed. Must always be given.
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[in] function
     * A \transient reference to the analytic function that is to be integrated.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the input function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::AnalyticFunctionIntegral helper class.
     */
    template<int max_der_, typename DataType_, typename Trafo_, typename Function_>
    typename AnalyticFunctionIntegral<DataType_, Function_>::Type integrate_analytic_function(
      DomainAssembler<Trafo_>& dom_asm, const Function_& function, const String& cubature)
    {
      DomainAssemblyAnalyticFunctionIntegralJob<DataType_, Function_, Trafo_, max_der_> job(function, dom_asm.get_trafo(), cubature);
      dom_asm.assemble(job);
      return job.result();
    }

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
    class DomainAssemblyDiscreteFunctionIntegralJob
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
        explicit Task(DomainAssemblyDiscreteFunctionIntegralJob& job) :
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
      explicit DomainAssemblyDiscreteFunctionIntegralJob(const Vector_& vector, const Space_& space, const String& cubature) :
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
    }; // class DomainAssemblyDiscreteFunctionIntegralJob<...>

    /**
     * \brief Assembles the integral of a discrete finite element function
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[in] vector
     * A \transient reference to the coefficient vector of the finite element function that is to be integrated.
     *
     * \param[in] space
     * A \transient reference to the finite element space that the coefficient vector belongs to.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the input function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::DiscreteFunctionIntegral helper class.
     */
    template<int max_der_, typename Vector_, typename Trafo_, typename Space_>
    typename DiscreteFunctionIntegral<Vector_, Space_>::Type integrate_discrete_function(
       DomainAssembler<Trafo_>& dom_asm, const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      DomainAssemblyDiscreteFunctionIntegralJob<Vector_, Space_, max_der_> job(vector, space, cubature);
      dom_asm.assemble(job);
      return job.result();
    }

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
    class DomainAssemblyErrorFunctionIntegralJob
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
        explicit Task(DomainAssemblyErrorFunctionIntegralJob& job) :
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
      explicit DomainAssemblyErrorFunctionIntegralJob(const Function_& function, const Vector_& vector, const Space_& space, const String& cubature) :
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
    }; // class DomainAssemblyErrorFunctionIntegralJob<...>

    /**
     * \brief Assembles the integral of an (analytic - discrete) error function
     *
     * \param[inout] dom_asm
     * A \transient reference to the domain assembler that is to be used for the assembly.
     *
     * \param[in] function
     * A \transient reference to the analytic function that is to be integrated.
     *
     * \param[in] vector
     * A \transient reference to the coefficient vector of the finite element function that is to be integrated.
     *
     * \param[in] space
     * A \transient reference to the finite element space that the coefficient vector belongs to.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the error function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::DiscreteFunctionIntegral helper class.
     */
    template<int max_der_, typename Function_, typename Vector_, typename Trafo_, typename Space_>
    typename DiscreteFunctionIntegral<Vector_, Space_>::Type integrate_error_function(
       DomainAssembler<Trafo_>& dom_asm, const Function_& function,
       const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(dom_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      DomainAssemblyErrorFunctionIntegralJob<Function_, Vector_, Space_, max_der_> job(function, vector, space, cubature);
      dom_asm.assemble(job);
      return job.result();
    }

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
    }; // class CellErrorFunctionIntegralJob<...>
  } // namespace Assembly
} // namespace FEAT
