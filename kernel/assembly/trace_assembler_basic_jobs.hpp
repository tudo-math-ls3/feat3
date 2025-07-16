// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/function_integral_info.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/intern/face_ref_trafo.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/congruency_trafo.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Basic assembly task base class for a single finite element space without pairwise assembly support
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Space_, TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags space_config_>
    class TraceAssemblyBasicTaskBase1
    {
    public:
      // we do not support pairwise assembly
      static constexpr bool assemble_pairwise = false;
      /// our assembly traits
      typedef Assembly::AsmTraits1<DataType_, Space_, trafo_config_, space_config_> AsmTraits;
      /// our data type
      typedef typename AsmTraits::DataType DataType;
      /// space type
      typedef typename AsmTraits::SpaceType SpaceType;
      /// trafo type
      typedef typename AsmTraits::TrafoType TrafoType;
      /// our shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// our facet type
      typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// the facet dimension; always equal to shape_dim-1
      static constexpr int facet_dim = FacetType::dimension;

      /// trafo evaluator type
      typedef typename AsmTraits::TrafoEvaluator TrafoEvaluator;
      /// trafo facet evaluator type
      typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;

      /// space evaluator type
      typedef typename AsmTraits::SpaceEvaluator SpaceEvaluator;

      /// include jacobian in facet trafo config (required for normal vector computation)
      static constexpr TrafoTags facet_trafo_config = facet_trafo_config_ | TrafoTags::jac_mat;

      /// trafo evaluation data
      typedef typename AsmTraits::TrafoEvalData TrafoEvalData;
      /// trafo facet evaluation data
      typedef typename TrafoFacetEvaluator::template ConfigTraits<facet_trafo_config>::EvalDataType TrafoFacetEvalData;
      /// space evaluation data
      typedef typename AsmTraits::SpaceEvalData SpaceEvalData;
      typedef typename AsmTraits::SpaceBasisData SpaceBasisData;

      /// cubature rule type
      typedef typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType CubatureRuleType;

    protected:
      /// the space
      const SpaceType& space;
      /// the trafo
      const TrafoType& trafo;
      /// the trafo evaluator
      TrafoEvaluator trafo_eval;
      /// the trafo facet evaluator
      TrafoFacetEvaluator trafo_facet_eval;
      /// the space evaluator
      SpaceEvaluator space_eval;
      /// the trafo evaluation data
      TrafoEvalData trafo_data;
      /// the trafo facet evaluation data
      TrafoFacetEvalData trafo_facet_data;
      /// the space evaluation data
      SpaceEvalData space_data;
      /// the space dof-mapping
      typename AsmTraits::DofMapping dof_mapping;

      /// the internal cell facet orientation code
      int cell_facet_ori;

      /// local facet trafo matrices and vectors
      Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
      Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
      Tiny::Vector<DataType, shape_dim> face_vec;
      Tiny::Vector<DataType, facet_dim> ori_vec;

      /// current cubature point on reference cell
      Tiny::Vector<DataType, shape_dim> cur_point;
      /// current cubature point on reference facet
      Tiny::Vector<DataType, facet_dim> cur_point_facet;

    public:
      explicit TraceAssemblyBasicTaskBase1(const SpaceType& space_) :
        space(space_),
        trafo(space.get_trafo()),
        trafo_eval(trafo),
        trafo_facet_eval(trafo),
        space_eval(space),
        trafo_data(),
        trafo_facet_data(),
        space_data(),
        dof_mapping(space),
        cell_facet_ori(0)
      {
        face_mat.format();
        ori_mat.format();
        face_vec.format();
        ori_vec.format();
      }

      virtual ~TraceAssemblyBasicTaskBase1() = default;

      void prepare(Index facet, Index cell, int local_facet, int facet_ori)
      {
        // compute facet trafo
        Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat, face_vec, local_facet);

        // compute orientation trafo
        Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat, ori_vec, facet_ori);

        // compute orientation of actual cell facet
        cell_facet_ori = Geometry::Intern::CongruencySampler<FacetType>::orientation(facet_ori)
          * Shape::ReferenceCell<ShapeType>::facet_orientation(local_facet);

        // prepare evaluators
        trafo_facet_eval.prepare(facet);
        trafo_eval.prepare(cell);
        space_eval.prepare(trafo_eval);

        // prepare dof mapping
        dof_mapping.prepare(cell);
      }

      void prepare_point(Tiny::Vector<DataType, facet_dim>& pt)
      {
        // set cubature point
        cur_point_facet = pt;

        // transform point to local facet on reference cell
        cur_point = (face_mat * ((ori_mat * cur_point_facet) + ori_vec)) + face_vec;

        // compute trafo data
        trafo_facet_eval(trafo_facet_data, cur_point_facet);
        trafo_eval(trafo_data, cur_point);

        // compute normal vector
        trafo_facet_data.normal = Tiny::orthogonal(trafo_facet_data.jac_mat).normalize();

        // ensure that the normal vector "points outside"
        if(cell_facet_ori < 0)
          trafo_facet_data.normal.negate();

        // copy normal to cell trafo data
        trafo_data.normal = trafo_facet_data.normal;

        // compute space data
        space_eval(space_data, trafo_data);
      }

      void finish()
      {
        // finish evaluator
        dof_mapping.finish();

        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();
        trafo_facet_eval.finish();
      }
    }; // class TraceAssemblyBasicTaskBase1<...>

    template<typename DataType_, typename TestSpace_, typename TrialSpace_,
      TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags test_config_, SpaceTags trial_config_>
    class TraceAssemblyBasicTaskBase2
    {
    public:
      // we do not support pairwise assembly
      static constexpr bool assemble_pairwise = false;
      /// our assembly traits
      typedef Assembly::AsmTraits2<DataType_, TestSpace_, TrialSpace_, trafo_config_, test_config_, trial_config_> AsmTraits;
      /// our data type
      typedef typename AsmTraits::DataType DataType;
      /// test space type
      typedef typename AsmTraits::TestSpaceType TestSpaceType;
      /// test space type
      typedef typename AsmTraits::TrialSpaceType TrialSpaceType;
      /// trafo type
      typedef typename AsmTraits::TrafoType TrafoType;
      /// our shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// our facet type
      typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// the facet dimension; always equal to shape_dim-1
      static constexpr int facet_dim = FacetType::dimension;

      /// trafo evaluator type
      typedef typename AsmTraits::TrafoEvaluator TrafoEvaluator;
      /// trafo facet evaluator type
      typedef typename TrafoType::template Evaluator<FacetType, DataType>::Type TrafoFacetEvaluator;

      /// test space evaluator type
      typedef typename AsmTraits::TestEvaluator TestEvaluator;
      /// trial space evaluator type
      typedef typename AsmTraits::TrialEvaluator TrialEvaluator;

      /// include jacobian in facet trafo config (required for normal vector computation)
      static constexpr TrafoTags facet_trafo_config = facet_trafo_config_ | TrafoTags::jac_mat;

      /// trafo evaluation data
      typedef typename AsmTraits::TrafoEvalData TrafoEvalData;
      /// trafo facet evaluation data
      typedef typename TrafoFacetEvaluator::template ConfigTraits<facet_trafo_config>::EvalDataType TrafoFacetEvalData;
      /// test evaluation data
      typedef typename AsmTraits::TestEvalData TestEvalData;
      typedef typename AsmTraits::TestBasisData TestBasisData;
      /// trial evaluation data
      typedef typename AsmTraits::TrialEvalData TrialEvalData;
      typedef typename AsmTraits::TrialBasisData TrialBasisData;

      /// cubature rule type
      typedef typename Intern::CubatureTraits<TrafoFacetEvaluator>::RuleType CubatureRuleType;

    protected:
      /// the test space
      const TestSpaceType& test_space;
      /// the trial space
      const TrialSpaceType& trial_space;
      /// the trafo
      const TrafoType& trafo;
      /// the trafo evaluator
      TrafoEvaluator trafo_eval;
      /// the trafo facet evaluator
      TrafoFacetEvaluator trafo_facet_eval;
      /// the test space evaluator
      TestEvaluator test_eval;
      /// the trial space evaluator
      TrialEvaluator trial_eval;
      /// the trafo evaluation data
      TrafoEvalData trafo_data;
      /// the trafo facet evaluation data
      TrafoFacetEvalData trafo_facet_data;
      /// the test space evaluation data
      TestEvalData test_data;
      /// the trial space evaluation data
      TrialEvalData trial_data;
      /// the space dof-mappings
      typename AsmTraits::TestDofMapping test_dof_mapping;
      typename AsmTraits::TrialDofMapping trial_dof_mapping;

      /// the internal cell facet orientation code
      int cell_facet_ori;

      /// local facet trafo matrices and vectors
      Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat;
      Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat;
      Tiny::Vector<DataType, shape_dim> face_vec;
      Tiny::Vector<DataType, facet_dim> ori_vec;

      /// current cubature point on reference cell
      Tiny::Vector<DataType, shape_dim> cur_point;
      /// current cubature point on reference facet
      Tiny::Vector<DataType, facet_dim> cur_point_facet;

    public:
      explicit TraceAssemblyBasicTaskBase2(const TestSpaceType& test_space_, const TrialSpaceType& trial_space_) :
        test_space(test_space_),
        trial_space(trial_space_),
        trafo(test_space.get_trafo()),
        trafo_eval(trafo),
        trafo_facet_eval(trafo),
        test_eval(test_space),
        trial_eval(trial_space),
        trafo_data(),
        trafo_facet_data(),
        test_data(),
        trial_data(),
        test_dof_mapping(test_space),
        trial_dof_mapping(trial_space),
        cell_facet_ori(0)
      {
        face_mat.format();
        ori_mat.format();
        face_vec.format();
        ori_vec.format();
      }

      virtual ~TraceAssemblyBasicTaskBase2() = default;

      void prepare(Index facet, Index cell, int local_facet, int facet_ori)
      {
        // compute facet trafo
        Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat, face_vec, local_facet);

        // compute orientation trafo
        Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat, ori_vec, facet_ori);

        // compute orientation of actual cell facet
        cell_facet_ori = Geometry::Intern::CongruencySampler<FacetType>::orientation(facet_ori)
          * Shape::ReferenceCell<ShapeType>::facet_orientation(local_facet);

        // prepare evaluators
        trafo_facet_eval.prepare(facet);
        trafo_eval.prepare(cell);
        test_eval.prepare(trafo_eval);
        trial_eval.prepare(trafo_eval);

        // prepare dof mapping
        test_dof_mapping.prepare(cell);
        trial_dof_mapping.prepare(cell);
      }

      void prepare_point(Tiny::Vector<DataType, facet_dim>& pt)
      {
        // set cubature point
        cur_point_facet = pt;

        // transform point to local facet on reference cell
        cur_point = (face_mat * ((ori_mat * cur_point_facet) + ori_vec)) + face_vec;

        // compute trafo data
        trafo_facet_eval(trafo_facet_data, cur_point_facet);
        trafo_eval(trafo_data, cur_point);

        // compute normal vector
        trafo_facet_data.normal = Tiny::orthogonal(trafo_facet_data.jac_mat).normalize();

        // ensure that the normal vector "points outside"
        if(cell_facet_ori < 0)
          trafo_facet_data.normal.negate();

        // copy normal to cell trafo data
        trafo_data.normal = trafo_facet_data.normal;

        // compute space data
        test_eval(test_data, trafo_data);
        trial_eval(trial_data, trafo_data);
      }

      void finish()
      {
        // finish evaluator
        trial_dof_mapping.finish();
        test_dof_mapping.finish();

        // finish evaluators
        trial_eval.finish();
        test_eval.finish();
        trafo_eval.finish();
        trafo_facet_eval.finish();
      }
    }; // class TraceAssemblyBasicTaskBase2<...>

    /**
     * \brief Basic Vector trace assembly task CRTP base-class
     *
     * This CRTP class can be used as a base-class for various types of vector assemblies tasks.
     * This class implements all the five interface member functions required by the TraceAssemblyJob
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
     * The cell trafo configuration tags required for the assembly
     *
     * \tparam facet_trafo_config_
     * The facet trafo configuration tags required for the assembly
     *
     * \tparam space_config_
     * The space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename Vector_, typename Space_,
      TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags space_config_>
    class TraceAssemblyBasicVectorTaskCRTP :
      public TraceAssemblyBasicTaskBase1<typename Vector_::DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_>
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      /// our base-class
      typedef TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_> BaseClass;
      /// our assembly traits
      typedef typename BaseClass::AsmTraits AsmTraits;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      using typename BaseClass::SpaceBasisData;

    protected:
      /// the vector that is to be assembled
      Vector_& vector;
      /// the cubature rule used for integration
      typename BaseClass::CubatureRuleType cubature_rule;
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
       * \param[in] cubature_factory_
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit TraceAssemblyBasicVectorTaskCRTP(Vector_& vector_, const Space_& space_,
        const Cubature::DynamicFactory& cubature_factory_, DataType alpha_) :
        BaseClass(space_),
        vector(vector_),
        cubature_rule(Cubature::ctor_factory, cubature_factory_),
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
       * \param[in] tau_f
       * The \transient facet trafo evaluation data for the current cubature point.
       *
       * \param[in] tau
       * The \transient trafo evaluation data for the current cubature point.
       */
      void set_point(TrafoFacetEvalData& tau_f, TrafoEvalData& tau);

      /**
       * \brief Evaluates the assembly for a basis function
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] psi
       * The \transient evaluation data of the current basis function.
       */
      void eval(ValueType& val, const DataType& weight, const SpaceBasisData& psi);
#endif // DOXYGEN

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local vector
        local_vector.format();

        // fetch number of local dofs
        const int num_loc_dofs = this->space_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // prepare point
          static_cast<Derived_&>(*this).prepare_point(cubature_rule.get_point(k));

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(this->trafo_facet_data, this->trafo_data);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // evaluate functional and integrate
            static_cast<Derived_&>(*this).eval(local_vector(i),
              this->trafo_facet_data.jac_det * cubature_rule.get_weight(k), this->space_data.phi[i]);
          } // continue with next test function
        } // continue with next cubature point
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        scatter_axpy(local_vector, this->dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class TraceAssemblyBasicVectorTaskCRTP<...>

    /**
     * \brief Vector trace assembly job for LinearFunctional implementations
     *
     * This class implements the TraceAssemblyJob interface to assemble a linearform that is
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
    class TraceAssemblyLinearFunctionalVectorJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = LinearFunctional_::trafo_config;
      static constexpr SpaceTags space_config = LinearFunctional_::test_config;

      class Task :
        public TraceAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, TrafoTags::none, space_config>
      {
      protected:
        /// our base-class typedef
        typedef TraceAssemblyBasicVectorTaskCRTP<Task, Vector_, Space_, trafo_config, TrafoTags::none, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename LinearFunctional_::template Evaluator<AsmTraits> func_eval;

      public:
        explicit Task(TraceAssemblyLinearFunctionalVectorJob& job) :
          BaseClass(job.vector, job.space, job.cubature_factory, job.alpha),
          func_eval(job.linear_functional)
        {
        }

        void prepare(Index facet, Index cell, int local_facet, int facet_ori)
        {
          BaseClass::prepare(facet, cell, local_facet, facet_ori);
          func_eval.prepare(this->trafo_eval);
        }

        void set_point(typename BaseClass::TrafoFacetEvalData& DOXY(tau_f), typename BaseClass::TrafoEvalData& tau)
        {
          func_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight, typename BaseClass::SpaceBasisData& psi)
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
      explicit TraceAssemblyLinearFunctionalVectorJob(const LinearFunctional_& linear_functional_,
        Vector_& vector_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        linear_functional(linear_functional_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class TraceAssemblyLinearFunctionalVectorJob<...>

    /**
     * \brief Assembles a linear functional into a vector
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
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
    void assemble_linear_functional_vector(TraceAssembler<Trafo_>& trace_asm,
      Vector_& vector, const LinFunc_& linear_functional, const Space_& space,
      const String& cubature, const typename Vector_::DataType alpha = typename Vector_::DataType(1))
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "trace assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector length");

      TraceAssemblyLinearFunctionalVectorJob<LinFunc_, Vector_, Space_> job(
        linear_functional, vector, space, cubature, alpha);
      trace_asm.assemble(job);
    }

    /**
     * \brief Basic Matrix trace assembly task CRTP base-class for identical test-/trial-space
     *
     * This CRTP class can be used as a base-class for various types of matrix assemblies tasks.
     * This class implements all the five interface member functions required by the TraceAssemblyJob
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
     * The finite element space to be used as test-space
     *
     * \tparam trafo_config_
     * The cell trafo configuration tags required for the assembly
     *
     * \tparam facet_trafo_config_
     * The facet trafo configuration tags required for the assembly
     *
     * \tparam space_config_
     * The space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename Matrix_, typename Space_,
      TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags space_config_>
    class TraceAssemblyBasicMatrixTaskCRTP1 :
      public TraceAssemblyBasicTaskBase1<typename Matrix_::DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_>
    {
    public:
      /// the data-type of the matrix
      typedef typename Matrix_::DataType DataType;
      /// the value-type of the matrix
      typedef typename Matrix_::ValueType ValueType;

      /// our base class
      typedef TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_> BaseClass;
      /// our assembly traits
      typedef typename BaseClass::AsmTraits AsmTraits;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      using typename BaseClass::SpaceBasisData;

    protected:
      /// the matrix that is to be assembled
      Matrix_& matrix;
      /// the cubature rule used for integration
      typename BaseClass::CubatureRuleType cubature_rule;
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
       * \param[in] cubature_factory_
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit TraceAssemblyBasicMatrixTaskCRTP1(Matrix_& matrix_, const Space_& space_,
        const Cubature::DynamicFactory& cubature_factory_, DataType alpha_) :
        BaseClass(space_),
        matrix(matrix_),
        cubature_rule(Cubature::ctor_factory, cubature_factory_),
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
       * The \transient cell trafo evaluation data for the current cubature point.
       *
       * \param[in] tau_f
       * The \transient facet trafo evaluation data for the current cubature point.
       */
      void set_point(TrafoFacetEvalData& tau_f, TrafoEvalData& tau);

      /**
       * \brief Evaluates the assembly for a test-/trial-basis function pair
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] psi
       * The \transient evaluation data of the current test basis function.
       *
       * \param[in] phi
       * The \transient evaluation data of the current trial basis function.
       */
      void eval(ValueType& val, const DataType& weight, const SpaceBasisData& psi, const SpaceBasisData& phi);
#endif // DOXYGEN

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local matrix
        local_matrix.format();

        // fetch number of local dofs
        const int num_loc_dofs = this->space_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // prepare point
          static_cast<Derived_&>(*this).prepare_point(cubature_rule.get_point(k));

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(this->trafo_facet_data, this->trafo_data);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_dofs; ++j)
            {
              // evaluate operator and integrate
              static_cast<Derived_&>(*this).eval(local_matrix(i,j),
                this->trafo_facet_data.jac_det * cubature_rule.get_weight(k),
                this->space_data.phi[j], this->space_data.phi[i]);
            } // continue with next trial function
          } // continue with next test function
        } // continue with next cubature point
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        scatter_axpy(local_matrix, this->dof_mapping, this->dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class TraceAssemblyBasicMatrixTaskCRTP1<...>

    /**
     * \brief Basic Matrix assembly task CRTP base-class for different test-/trial-space
     *
     * This CRTP class can be used as a base-class for various types of matrix assemblies tasks.
     * This class implements all the five interface member functions required by the TraceAssemblyJob
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
     * The cell trafo configuration tags required for the assembly
     *
     * \tparam facet_trafo_config_
     * The facet trafo configuration tags required for the assembly
     *
     * \tparam test_config_
     * The test-space configuration tags required for the assembly
     *
     * \tparam trial_config_
     * The trial-space configuration tags required for the assembly
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename Matrix_, typename TestSpace_, typename TrialSpace_,
      TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags test_config_, SpaceTags trial_config_>
    class TraceAssemblyBasicMatrixTaskCRTP2 :
      public TraceAssemblyBasicTaskBase2<typename Matrix_::DataType, TestSpace_, TrialSpace_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, test_config_, trial_config_>
    {
    public:
      /// the data-type of the matrix
      typedef typename Matrix_::DataType DataType;
      /// the value-type of the matrix
      typedef typename Matrix_::ValueType ValueType;

      /// our base class
      typedef TraceAssemblyBasicTaskBase2<typename Matrix_::DataType, TestSpace_, TrialSpace_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, test_config_, trial_config_> BaseClass;
      /// our assembly traits
      typedef typename BaseClass::AsmTraits AsmTraits;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      using typename BaseClass::TestBasisData;
      using typename BaseClass::TrialBasisData;

    protected:
      /// the matrix that is to be assembled
      Matrix_& matrix;
      /// the cubature rule used for integration
      typename BaseClass::CubatureRuleType cubature_rule;
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
       * A \resident reference to the finite element test space to be used for the assembly
       *
       * \param[in] trial_space_
       * A \resident reference to the finite element trial space to be used for the assembly
       *
       * \param[in] cubature_factory_
       * A \transient reference to the cubature factory that is used for the cubature rule.
       *
       * \param[in] alpha_
       * A scaling factor for the assembly.
       */
      explicit TraceAssemblyBasicMatrixTaskCRTP2(Matrix_& matrix_, const TestSpace_& test_space_,
        const TestSpace_& trial_space_, const Cubature::DynamicFactory& cubature_factory_, DataType alpha_) :
        BaseClass(test_space_, trial_space_),
        matrix(matrix_),
        cubature_rule(Cubature::ctor_factory, cubature_factory_),
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
       * \param[in] tau_f
       * The \transient facet trafo evaluation data for the current cubature point.
       *
       * \param[in] tau
       * The \transient cell trafo evaluation data for the current cubature point.
       */
      void set_point(TrafoFacetEvalData& tau_f, TrafoEvalData& tau);

      /**
       * \brief Evaluates the assembly for a test-/trial-basis function pair
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
      void eval(ValueType& val, const DataType& weight, const TrialBasisData& phi, const TestBasisData& psi),
#endif // DOXYGEN

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local matrix
        local_matrix.format();

        // fetch number of local dofs
        const int num_loc_test_dofs = this->test_eval.get_num_local_dofs();
        const int num_loc_trial_dofs = this->trial_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // prepare point
          static_cast<Derived_&>(*this).prepare_point(cubature_rule.get_point(k));

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(this->trafo_facet_data, this->trafo_data);

          // test function loop
          for(int i(0); i < num_loc_test_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_trial_dofs; ++j)
            {
              // evaluate operator and integrate
              static_cast<Derived_&>(*this).eval(local_matrix(i,j),
                this->trafo_facet_data.jac_det * cubature_rule.get_weight(k),
                this->trial_data.phi[j], this->test_data.phi[i]);
            } // continue with next trial function
          } // continue with next test function
        } // continue with next cubature point
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        scatter_axpy(local_matrix, this->test_dof_mapping, this->trial_dof_mapping, scatter_alpha);
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class TraceAssemblyBasicMatrixTaskCRTP2<...>

    /**
     * \brief Matrix assembly job for BilinearOperator implementations and identical test-/trial-spaces
     *
     * This class implements the TraceAssemblyJob interface to assemble a bilinearform that is
     * provided as a BilinearOperator implementation.
     *
     * \tparam BilinearOperator_
     * The bilinear operator that is to be assembled
     *
     * \tparam Matrix_
     * The vector that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename BilinearOperator_, typename Matrix_, typename Space_>
    class TraceAssemblyBilinearOperatorMatrixJob1
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags space_config = BilinearOperator_::test_config | BilinearOperator_::trial_config;

      class Task :
        public TraceAssemblyBasicMatrixTaskCRTP1<Task, Matrix_, Space_, trafo_config, TrafoTags::none, space_config>
      {
      protected:
        /// our base-class typedef
        typedef TraceAssemblyBasicMatrixTaskCRTP1<Task, Matrix_, Space_, trafo_config, TrafoTags::none, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(TraceAssemblyBilinearOperatorMatrixJob1& job) :
          BaseClass(job.matrix, job.space, job.cubature_factory, job.alpha),
          oper_eval(job.bilinear_operator)
        {
        }

        void prepare(Index facet, Index cell, int local_facet, int facet_ori)
        {
          BaseClass::prepare(facet, cell, local_facet, facet_ori);
          oper_eval.prepare(this->trafo_eval);
        }

        void set_point(typename BaseClass::TrafoFacetEvalData& DOXY(tau_f), typename BaseClass::TrafoEvalData& tau)
        {
          oper_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight, typename AsmTraits::TrialBasisData& phi, typename AsmTraits::TestBasisData& psi)
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
      explicit TraceAssemblyBilinearOperatorMatrixJob1(const BilinearOperator_& bilinear_operator_,
        Matrix_& matrix_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        bilinear_operator(bilinear_operator_),
        matrix(matrix_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class TraceAssemblyBilinearOperatorMatrixJob1<...>

    /**
     * \brief Assembles a bilinear operator into a matrix with identical test- and trial-spaces.
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
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
    void assemble_bilinear_operator_matrix_1(TraceAssembler<Trafo_>& trace_asm, Matrix_& matrix,
      const BilOp_& bilinear_operator, const Space_& space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "trace assembler and space have different trafos");
      XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix row count");

      TraceAssemblyBilinearOperatorMatrixJob1<BilOp_, Matrix_, Space_> job(
        bilinear_operator, matrix, space, cubature, alpha);
      trace_asm.assemble(job);
    }

    /**
     * \brief Matrix assembly job for BilinearOperator implementations and different test-/trial-spaces
     *
     * This class implements the TraceAssemblyJob interface to assemble a bilinearform that is
     * provided as a BilinearOperator implementation.
     *
     * \tparam BilinearOperator_
     * The bilinear operator that is to be assembled
     *
     * \tparam Matrix_
     * The vector that is to be assembled.
     *
     * \tparam TestSpace_
     * The finite element space to be used as test-space
     *
     * \tparam TrialSpace_
     * The finite element space to be used as trial-space
     *
     * \author Peter Zajac
     */
    template<typename BilinearOperator_, typename Matrix_, typename TestSpace_, typename TrialSpace_>
    class TraceAssemblyBilinearOperatorMatrixJob2
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = BilinearOperator_::trafo_config;
      static constexpr SpaceTags test_config  = BilinearOperator_::test_config;
      static constexpr SpaceTags trial_config = BilinearOperator_::trial_config;

      class Task :
        public TraceAssemblyBasicMatrixTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, TrafoTags::none, test_config, trial_config>
      {
      protected:
        /// our base-class typedef
        typedef TraceAssemblyBasicMatrixTaskCRTP2<Task, Matrix_, TestSpace_, TrialSpace_, trafo_config, TrafoTags::none, test_config, trial_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// the bilinear operator evaluator
        typename BilinearOperator_::template Evaluator<AsmTraits> oper_eval;

      public:
        explicit Task(TraceAssemblyBilinearOperatorMatrixJob2& job) :
          BaseClass(job.matrix, job.test_space, job.trial_space, job.cubature_factory, job.alpha),
          oper_eval(job.bilinear_operator)
        {
        }

        void prepare(Index facet, Index cell, int local_facet, int facet_ori)
        {
          BaseClass::prepare(facet, cell, local_facet, facet_ori);
          oper_eval.prepare(this->trafo_eval);
        }

        void set_point(typename BaseClass::TrafoFacetEvalData& DOXY(tau_f), typename BaseClass::TrafoEvalData& tau)
        {
          oper_eval.set_point(tau);
        }

        void eval(ValueType& val, const DataType& weight, typename AsmTraits::TrialBasisData& phi, typename AsmTraits::TestBasisData& psi)
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
      explicit TraceAssemblyBilinearOperatorMatrixJob2(const BilinearOperator_& bilinear_operator_,
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
    }; // class TraceAssemblyBilinearOperatorMatrixJob2<...>

    /**
     * \brief Assembles a bilinear operator into a matrix with different test- and trial-spaces.
     *
     * \param[inout] trace_asm
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
    void assemble_bilinear_operator_matrix_2(TraceAssembler<Trafo_>& trace_asm,
      Matrix_& matrix, const BilOp_& bilinear_operator, const TestSpace_& test_space,
      const TrialSpace_& trial_space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(trace_asm.get_trafo() == test_space.get_trafo(), "trace assembler and test space have different trafos");
      XASSERTM(trace_asm.get_trafo() == trial_space.get_trafo(), "trace assembler and trial space have different trafos");
      XASSERTM(matrix.columns() == trial_space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == test_space.get_num_dofs(), "invalid matrix row count");

      TraceAssemblyBilinearOperatorMatrixJob2<BilOp_, Matrix_, TestSpace_, TrialSpace_> job(
        bilinear_operator, matrix, test_space, trial_space, cubature, alpha);
      trace_asm.assemble(job);
    }

    /**
     * \brief Assembly job for the trace integration of an analytic function
     *
     * \tparam DataType_
     * The (scalar) data type in which the assembly is to be performed. Must always be given.
     *
     * \tparam Function_
     * The type of the analytic function u that is to be integrated.
     *
     * \tparam Trafo_
     * The transformation that discretizes the integration domain.
     *
     * \tparam max_der_
     * The maximum derivative of the function that is to be integrated.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Function_, typename Trafo_, int max_der_>
    class TraceAssemblyAnalyticFunctionIntegralJob
    {
    public:
      typedef DataType_ DataType;

      /// our function integral type
      typedef typename AnalyticFunctionIntegral<DataType, Function_>::Type FunctionIntegralType;

      class Task
      {
      public:
        /// our base-class typedef
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;

        typedef typename Trafo_::ShapeType ShapeType;
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::ShapeType FacetType;

        /// trafo evaluator
        typedef typename Trafo_::template Evaluator<FacetType, DataType>::Type TrafoEvaluator;
        typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;

        /// declare our analytic eval traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = true;

      protected:
        /// trafo evaluator
        TrafoEvaluator trafo_eval;
        /// trafo eval data type
        TrafoEvalData trafo_data;
        /// the cubature rule
        typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType cubature_rule;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(TraceAssemblyAnalyticFunctionIntegralJob& job) :
          trafo_eval(job.trafo),
          trafo_data(),
          cubature_rule(Cubature::ctor_factory, job.cubature_factory),
          func_eval(job.function),
          loc_integral(),
          job_integral(job.integral)
        {
        }

        void prepare(Index facet, Index DOXY(cell), int DOXY(local_facet), int DOXY(facet_ori))
        {
          trafo_eval.prepare(facet);
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
      const Trafo_& trafo;
      /// the function to be integrated
      const Function_& function;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the function integral
      FunctionIntegralType integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function_
       * A \resident reference to the analytic function that is to be integrated.
       *
       * \param[in] trafo_
       * A \resident reference to the trafo representing the domain over which to integrate.
       *
       * \param[in] cubature_
       * The name of the cubature rule to use for integration.
       */
      explicit TraceAssemblyAnalyticFunctionIntegralJob(const Function_& function_, const Trafo_& trafo_, const String& cubature_) :
        trafo(trafo_),
        function(function_),
        cubature_factory(cubature_),
        integral()
      {
        integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      FunctionIntegralType& result()
      {
        return integral;
      }
    }; // class TraceAssemblyAnalyticFunctionIntegralJob<...>

    /**
     * \brief Assembles the trace integral of an analytic function
     *
     * \tparam max_der_
     * The maximum derivative of the function that is to be integrated.
     *
     * \tparam DataType_
     * The (scalar) data type in which the assembly is to be performed. Must always be given.
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
     *
     * \param[in] function
     * A \transient reference to the analytic function that is to be integrated.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \returns
     * The integrals of the analytic function stored in an Assembly::FunctionIntegralInfo instance,
     * whose exact type is determined by the Assembly::AnalyticFunctionIntegral helper class.
     */
    template<int max_der_, typename DataType_, typename Function_, typename Trafo_>
    typename AnalyticFunctionIntegral<DataType_, Function_>::Type integrate_analytic_function(
      TraceAssembler<Trafo_>& trace_asm, const Function_& function, const String& cubature)
    {
      TraceAssemblyAnalyticFunctionIntegralJob<DataType_, Function_, Trafo_, max_der_> job(function, trace_asm.get_trafo(), cubature);
      trace_asm.assemble(job);
      return job.result();
    }

    /**
     * \brief Assembly job for the trace integration of a discrete finite element function
     *
     * \tparam Vector_
     * The vector that represents the discrete function to be integrated.
     *
     * \tparam Space_
     * The finite element space that is to be used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative of the function that is to be integrated.
     *
     * \author Peter Zajac
     */
    template<typename Vector_, typename Space_, int max_der_>
    class TraceAssemblyDiscreteFunctionIntegralJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::none;
      static constexpr SpaceTags space_config = SpaceTags::value |
        (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
        (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);
      static constexpr TrafoTags facet_trafo_config = TrafoTags::jac_det;

      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

      class Task :
        public TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config>
      {
      public:
        /// our base-class typedef
        typedef TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = true;

        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::TrafoFacetEvalData;

      protected:
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the cubature rule used for integration
        typename BaseClass::CubatureRuleType cubature_rule;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector scatter object
        typename Vector_::GatherAxpy gather_axpy;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(TraceAssemblyDiscreteFunctionIntegralJob& job) :
          BaseClass(job.space),
          vector(job.vector),
          cubature_rule(Cubature::ctor_factory, job.cubature_factory),
          local_vector(),
          gather_axpy(vector),
          loc_integral(),
          job_integral(job.integral)
        {
        }

        /// \brief Performs the local assembly.
        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, this->dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = this->space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // prepare point
            this->prepare_point(cubature_rule.get_point(k));

            // do the dirty work
            Intern::DiscFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * this->trafo_facet_data.jac_det, this->space_data, local_vector, num_loc_dofs);
          } // continue with next cubature point
        }

        /// \brief Scatters the local assembly.
        void scatter()
        {
          // nothing to do here
        }

        /// \brief Finalizes the assembly.
        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// a reference to the vector that is to be integrated
      const Vector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the function integral
      FunctionIntegralType integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector_
       * A \resident reference to the coefficient vector of the discrete function.
       *
       * \param[in] space_
       * A \resident reference to the finite element space used for the discretization
       *
       * \param[in] cubature_
       * The name of the cubature rule to use for integration.
       */
      explicit TraceAssemblyDiscreteFunctionIntegralJob(const Vector_& vector_, const Space_& space_, String cubature_):
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        integral()
      {
        integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      FunctionIntegralType& result()
      {
        return integral;
      }
    }; // class TraceAssemblyDiscreteFunctionIntegralJob<...>

    /**
     * \brief Assembles the integral of a discrete finite element function
     *
     * \tparam max_der_
     * The maximum derivative of the function that is to be integrated.
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
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
      TraceAssembler<Trafo_>& trace_asm, const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "domain assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      TraceAssemblyDiscreteFunctionIntegralJob<Vector_, Space_, max_der_> job(vector, space, cubature);
      trace_asm.assemble(job);
      return job.result();
    }

    /**
     * \brief Assembly job for the trace integration of a analytic vs discrete error function
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
     * The finite element space to be used for the discretization
     *
     * \tparam max_der_
     * The maximum derivative of the function that is to be integrated.
     *
     * \author Peter Zajac
     */
    template<typename Function_, typename Vector_, typename Space_, int max_der_>
    class TraceAssemblyErrorFunctionIntegralJob
    {
    public:
      typedef typename Vector_::DataType DataType;
      typedef typename Vector_::ValueType ValueType;

      /// make sure that function and vector are compatible
      static_assert(Intern::ErrCompatHelper<Function_, Vector_>::valid, "function and vector are incompatible");

      /// our function integral type
      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

      static constexpr TrafoTags trafo_config = TrafoTags::none;
      static constexpr SpaceTags space_config = SpaceTags::value |
        (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
        (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);
      static constexpr TrafoTags facet_trafo_config = TrafoTags::img_point | TrafoTags::jac_det;

      class Task :
        public TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config>
      {
      public:
        /// our base-class typedef
        typedef TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        /// declare our analytic eval traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = true;

        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::TrafoFacetEvalData;

      protected:
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the cubature rule used for integration
        typename BaseClass::CubatureRuleType cubature_rule;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector scatter object
        typename Vector_::GatherAxpy gather_axpy;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(TraceAssemblyErrorFunctionIntegralJob& job) :
          BaseClass(job.space),
          func_eval(job.function),
          vector(job.vector),
          cubature_rule(Cubature::ctor_factory, job.cubature_factory),
          local_vector(),
          gather_axpy(vector),
          loc_integral(),
          job_integral(job.integral)
        {
        }

        /// \brief Performs the local assembly.
        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, this->dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = this->space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // prepare point
            this->prepare_point(cubature_rule.get_point(k));

            // do the dirty work
            Intern::ErrFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * this->trafo_facet_data.jac_det, func_eval,
              this->trafo_facet_data.img_point, this->space_data, local_vector, num_loc_dofs);

          } // continue with next cubature point
        }

        /// \brief Scatters the local assembly.
        void scatter()
        {
          // nothing to do here
        }

        /// \brief Finalizes the assembly.
        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the function to be integrated
      const Function_& function;
      /// a reference to the vector that is to be integrated
      const Vector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the function integral
      FunctionIntegralType integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function_
       * A \resident reference to the analytic function
       *
       * \param[in] vector_
       * A \resident reference to the coefficient vector of the discrete function.
       *
       * \param[in] space_
       * A \resident reference to the finite element space used for the discretization
       *
       * \param[in] cubature_
       * The name of the cubature rule to use for integration.
       */
      explicit TraceAssemblyErrorFunctionIntegralJob(const Function_& function_, const Vector_& vector_, const Space_& space_, String cubature_):
        function(function_),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_),
        integral()
      {
        integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      FunctionIntegralType& result()
      {
        return integral;
      }
    }; // class TraceAssemblyErrorFunctionIntegralJob<...>

    /**
     * \brief Assembles the trace integral of an (analytic - discrete) error function
     *
     * \tparam max_der_
     * The maximum derivative of the error function that is to be integrated.
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
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
      TraceAssembler<Trafo_>& trace_asm, const Function_& function,
      const Vector_& vector, const Space_& space, const String& cubature)
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "trace assembler and space have different trafos");
      XASSERTM(vector.size() == space.get_num_dofs(), "invalid coefficient vector length");
      TraceAssemblyErrorFunctionIntegralJob<Function_, Vector_, Space_, max_der_> job(function, vector, space, cubature);
      trace_asm.assemble(job);
      return job.result();
    }

    /**
     * \brief Basic Stokes Vector analysis task CRTP base-class
     *
     * This class can be used as a base class for some sort of post-processing assembly of a Stokes solution vector.
     * The most prominent use case of this class is the computation of body forces on an obstacle.
     *
     * \tparam Derived_
     * The most derived CRTP class
     *
     * \tparam VectorVelo_
     * The type of the velocity vector that is to be assembled
     *
     * \tparam VectorPres_
     * The type of the pressure vector that is to be assembled
     *
     * \tparam SpaceVelo_
     * The finite element velocity space
     *
     * \tparam SpacePres_
     * The finite element pressure space
     *
     * \tparam trafo_config_
     * The cell trafo configuration tags required for the assembly
     *
     * \tparam facet_trafo_config_
     * The facet trafo configuration tags required for the assembly
     *
     * \tparam space_velo_config_
     * The velocity space configuration tags required for the assembly
     *
     * \tparam space_pres_config_
     * The pressure space configuration tags required for the assembly
     *
     * \note
     * This class (ab)uses the base class for different trial/test functions to define the assembly traits and all the
     * required helper structures for the velocity/pressure finite element space pair.
     *
     * \author Peter Zajac
     */
    template<typename Derived_, typename VectorVelo_, typename VectorPres_, typename SpaceVelo_, typename SpacePres_,
      TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags space_velo_config_, SpaceTags space_pres_config_>
    class TraceAssemblyStokesVectorAnalysisTaskCRTP :
      public TraceAssemblyBasicTaskBase2<typename VectorVelo_::DataType, SpaceVelo_, SpacePres_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_velo_config_, space_pres_config_>
    {
    public:
      /// the data-type of the vector
      typedef typename VectorVelo_::DataType DataType;

      /// our base-class
      typedef TraceAssemblyBasicTaskBase2<DataType, SpaceVelo_, SpacePres_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_velo_config_, space_pres_config_> BaseClass;
      /// our assembly traits
      typedef typename BaseClass::AsmTraits AsmTraits;

      /// the value-type of the vectors
      typedef typename VectorVelo_::ValueType VeloValueType;
      typedef typename VectorPres_::ValueType PresValueType;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      typedef typename TrafoEvalData::EvalTraits TrafoEvalTraits;

      static constexpr int max_local_velo_dofs = AsmTraits::max_local_test_dofs;
      static constexpr int max_local_pres_dofs = AsmTraits::max_local_trial_dofs;

      typedef Space::BasisData<Space::StandardVectorEvalTraits<TrafoEvalTraits, max_local_velo_dofs, DataType>, space_velo_config_> VeloData;
      typedef Space::BasisData<Space::StandardScalarEvalTraits<TrafoEvalTraits, max_local_pres_dofs, DataType>, space_pres_config_> PresData;

    protected:
      /// the velocity vector that is to be analyzed
      const VectorVelo_& vector_velo;
      /// the pressure vector that is to be analyzed
      const VectorPres_& vector_pres;
      /// the cubature rule used for integration
      typename BaseClass::CubatureRuleType cubature_rule;
      /// the local velocity vector to be assembled
      typename AsmTraits::template TLocalTestVector<VeloValueType> local_vector_velo;
      /// the local pressure vector to be assembled
      typename AsmTraits::template TLocalTrialVector<PresValueType> local_vector_pres;
      /// the velocity vector object
      typename VectorVelo_::GatherAxpy gather_axpy_velo;
      /// the pressure vector object
      typename VectorPres_::GatherAxpy gather_axpy_pres;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector_velo_, vector_pres_
       * Two \resident references to the velocity and pressure vectors that are to be analyzed
       *
       * \param[in] space_velo_, space_pres_
       * Two \resident references to the velocity and pressure finite element spaces
       *
       * \param[in] cubature_factory_
       * A \transient reference to the cubature factory that is used for the cubature rule.
       */
      explicit TraceAssemblyStokesVectorAnalysisTaskCRTP(
        const VectorVelo_& vector_velo_, const VectorPres_& vector_pres_,
        const SpaceVelo_& space_velo_, const SpacePres_& space_pres_,
        const Cubature::DynamicFactory& cubature_factory_) :
        BaseClass(space_velo_, space_pres_),
        vector_velo(vector_velo_),
        vector_pres(vector_pres_),
        cubature_rule(Cubature::ctor_factory, cubature_factory_),
        gather_axpy_velo(vector_velo),
        gather_axpy_pres(vector_pres)
      {
      }

#ifdef DOXYGEN
      /**
       * \brief Sets the current cubature point
       *
       * \attention
       * This function must be implemented in the derived CRTP class!
       *
       * \param[in] tau_f
       * The \transient facet trafo evaluation data for the current cubature point.
       *
       * \param[in] tau
       * The \transient trafo evaluation data for the current cubature point.
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] velo
       * The \transient evaluation data of the velocity vector field in the current cubature point.
       *
       * \param[in] pres
       * The \transient evaluation data of the pressure function in the current cubature point.
       */
      void eval(TrafoFacetEvalData& tau_f, TrafoEvalData& tau, const DataType& weight,
        const VeloData& velo, const PresData& pres);
#endif // DOXYGEN

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // gather local vectors
        local_vector_velo.format();
        local_vector_pres.format();
        gather_axpy_velo(local_vector_velo, this->test_dof_mapping);
        gather_axpy_pres(local_vector_pres, this->trial_dof_mapping);

        // fetch number of local dofs
        const int num_loc_dofs_velo = this->test_eval.get_num_local_dofs();
        const int num_loc_dofs_pres = this->trial_eval.get_num_local_dofs();

        VeloData velo_data;
        PresData pres_data;

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // prepare point
          static_cast<Derived_&>(*this).prepare_point(cubature_rule.get_point(k));

          velo_data.format();
          pres_data.format();

          // compute velocity data
          for(int i(0); i < num_loc_dofs_velo; ++i)
          {
            if constexpr (*(space_velo_config_ & SpaceTags::value))
              velo_data.value.axpy(this->test_data.phi[i].value, local_vector_velo(i));
            if constexpr (*(space_velo_config_ & SpaceTags::ref_value))
              velo_data.ref_value.axpy(this->test_data.phi[i].ref_value, local_vector_velo(i));
            if constexpr (*(space_velo_config_ & SpaceTags::grad))
              velo_data.grad.add_outer_product(local_vector_velo(i), this->test_data.phi[i].grad);
            if constexpr (*(space_velo_config_ & SpaceTags::ref_grad))
              velo_data.ref_grad.add_outer_product(local_vector_velo(i), this->test_data.phi[i].ref_grad);
            if constexpr (*(space_velo_config_ & SpaceTags::hess))
              velo_data.hess.add_vec_mat_outer_product(local_vector_velo(i), this->test_data.phi[i].hess);
            if constexpr (*(space_velo_config_ & SpaceTags::ref_hess))
              velo_data.ref_hess.add_vec_mat_outer_product(local_vector_velo(i), this->test_data.phi[i].ref_hess);
          }

          // compute velocity data
          for(int i(0); i < num_loc_dofs_pres; ++i)
          {
            if constexpr (*(space_pres_config_ & SpaceTags::value))
              pres_data.value += local_vector_pres(i) * this->trial_data.phi[i].value;
            if constexpr (*(space_pres_config_ & SpaceTags::ref_value))
              pres_data.ref_value += local_vector_pres(i) * this->trial_data.phi[i].ref_value;
            if constexpr (*(space_pres_config_ & SpaceTags::grad))
              pres_data.grad.axpy(local_vector_pres(i), this->trial_data.phi[i].grad);
            if constexpr (*(space_pres_config_ & SpaceTags::ref_grad))
              pres_data.ref_grad.axpy(local_vector_pres(i), this->trial_data.phi[i].ref_grad);
            if constexpr (*(space_pres_config_ & SpaceTags::hess))
              pres_data.hess.axpy(local_vector_pres(i), this->trial_data.phi[i].hess);
            if constexpr (*(space_pres_config_ & SpaceTags::ref_hess))
              pres_data.ref_hess.axpy(local_vector_pres(i), this->trial_data.phi[i].ref_hess);
          }

          // evaluate
          static_cast<Derived_&>(*this).eval(this->trafo_facet_data, this->trafo_data,
            cubature_rule.get_weight(k) * this->trafo_facet_data.jac_det, velo_data, pres_data);
        } // continue with next cubature point
      }
    }; // class TraceAssemblyStokesVectorAnalysisTaskCRTP<...>

    /**
     * \brief Assembly job for the body forces computation of a Stokes solution vector
     *
     * \tparam VectorVelo_, VectorPres_
     * The vector types of the velocity and pressure components, respectively.
     *
     * \tparam SpaceVelo_, SpacePres_
     * The finite element spaces for the velocity and pressure discretizations
     *
     * This class implements the body force computation as proposed by Giles et al in \cite GLLS96
     *
     * \author Peter Zajac
     */
    template<typename VectorVelo_, typename VectorPres_, typename SpaceVelo_, typename SpacePres_>
    class TraceAssemblyStokesBodyForceAssemblyJob
    {
    public:
      typedef typename VectorVelo_::DataType DataType;
      typedef typename VectorVelo_::ValueType VeloValueType;
      typedef typename VectorPres_::ValueType PresValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::none;
      static constexpr SpaceTags space_velo_config = SpaceTags::value | SpaceTags::grad;
      static constexpr SpaceTags space_pres_config = SpaceTags::value;
      static constexpr TrafoTags facet_trafo_config = TrafoTags::jac_det | TrafoTags::normal;

      static constexpr int dim = SpaceVelo_::shape_dim;

      typedef Tiny::Matrix<DataType, 3, 2> RawForces;

      class Task :
        public TraceAssemblyStokesVectorAnalysisTaskCRTP<Task, VectorVelo_, VectorPres_, SpaceVelo_, SpacePres_, trafo_config, facet_trafo_config, space_velo_config, space_pres_config>
      {
      public:
        /// our base-class typedef
        typedef TraceAssemblyStokesVectorAnalysisTaskCRTP<Task, VectorVelo_, VectorPres_, SpaceVelo_, SpacePres_, trafo_config, facet_trafo_config, space_velo_config, space_pres_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = true;

        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::TrafoFacetEvalData;
        using typename BaseClass::VeloData;
        using typename BaseClass::PresData;

      protected:
        RawForces& job_raw_forces;
        RawForces raw_forces;

      public:
        explicit Task(TraceAssemblyStokesBodyForceAssemblyJob& job) :
          BaseClass(job.vector_velo, job.vector_pres, job.space_velo, job.space_pres, job.cubature_factory),
          job_raw_forces(job.raw_forces),
          raw_forces()
        {
          raw_forces.format();
        }

        void eval(TrafoFacetEvalData& tau_f, TrafoEvalData& DOXY(tau), DataType weight, const VeloData& velo, const PresData& pres)
        {
          this->_eval(weight, tau_f.normal, velo.grad, pres.value);
        }

        /// \brief Scatters the local assembly.
        void scatter()
        {
          // nothing to do here
        }

        /// \brief Finalizes the assembly.
        void combine()
        {
          job_raw_forces += raw_forces;
        }

      protected:
        /// 2D version
        void _eval(const DataType omega, const Tiny::Vector<DataType, 2, 2>& n,
          const Tiny::Matrix<DataType, 2, 2, 2, 2>& grad_v, const DataType val_p)
        {
          raw_forces(0, 0) -= omega * (DataType(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1]);
          raw_forces(1, 0) -= omega * (DataType(2) * grad_v(1,1) * n[1] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
          raw_forces(0, 1) -= omega * val_p * n[0];
          raw_forces(1, 1) -= omega * val_p * n[1];
        }

        /// 3D version
        void _eval(const DataType omega, const Tiny::Vector<DataType, 3, 3>& n,
          const Tiny::Matrix<DataType, 3, 3, 3, 3>& grad_v, const DataType val_p)
        {
          raw_forces(0, 0) -= omega * (DataType(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1] + (grad_v(0, 2) + grad_v(2, 0)) * n[2]);
          raw_forces(1, 0) -= omega * (DataType(2) * grad_v(1,1) * n[1] + (grad_v(1, 2) + grad_v(2, 1)) * n[2] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
          raw_forces(2, 0) -= omega * (DataType(2) * grad_v(2,2) * n[2] + (grad_v(2, 0) + grad_v(0, 2)) * n[0] + (grad_v(2, 1) + grad_v(1, 2)) * n[1]);
          raw_forces(0, 1) -= omega * val_p * n[0];
          raw_forces(1, 1) -= omega * val_p * n[1];
          raw_forces(2, 1) -= omega * val_p * n[2];
        }
      }; // class Task

    protected:
      const VectorVelo_& vector_velo;
      const VectorPres_& vector_pres;
      const SpaceVelo_& space_velo;
      const SpacePres_& space_pres;
      Cubature::DynamicFactory cubature_factory;
      RawForces raw_forces;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector_velo_
       * A \resident reference to the velocity coefficient vector
       *
       * \param[in] vector_pres_
       * A \resident reference to the pressure coefficient vector
       *
       * \param[in] space_velo_
       * A \resident reference to the finite element velocity space
       *
       * \param[in] space_pres_
       * A \resident reference to the finite element pressure space
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit TraceAssemblyStokesBodyForceAssemblyJob(const VectorVelo_& vector_velo_, const VectorPres_& vector_pres_,
        const SpaceVelo_& space_velo_, const SpacePres_& space_pres_, String cubature_):
        vector_velo(vector_velo_),
        vector_pres(vector_pres_),
        space_velo(space_velo_),
        space_pres(space_pres_),
        cubature_factory(cubature_),
        raw_forces()
      {
        raw_forces.format();
      }

      /// Returns the raw drag force coefficient for a given viscosity parameter
      DataType drag(DataType nu) const
      {
        return nu * raw_forces(0, 0) - raw_forces(0, 1);
      }

      /// Returns the raw lift  force coefficient for a given viscosity parameter
      DataType lift(DataType nu) const
      {
        return nu * raw_forces(1, 0) - raw_forces(1, 1);
      }

      /// Returns the raw side force coefficient for a given viscosity parameter
      DataType side(DataType nu) const
      {
        return nu * raw_forces(2, 0) - raw_forces(2, 1);
      }

      /// Synchronizes the forces over a communicator
      void sync(const Dist::Comm& comm)
      {
        comm.allreduce(&raw_forces.v[0][0], &raw_forces.v[0][0], std::size_t(6), Dist::op_sum);
      }
    }; // class TraceAssemblyStokesBodyForceAssemblyJob<...>
  } // namespace Assembly
} // namespace FEAT
