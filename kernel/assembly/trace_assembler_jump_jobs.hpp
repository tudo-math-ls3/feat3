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
    /// \cond internal
    namespace Intern
    {
      template<int max_>
      class CommonDofMap
      {
      public:
        Index glob_dofs[max_];
        int loc1_dofs[max_];
        int loc2_dofs[max_];
        int num_dofs;

        explicit CommonDofMap() : num_dofs(0) {}

        void clear()
        {
          num_dofs = 0;
        }

        int push_1(Index dof, int loc)
        {
          for(int i(0); i < num_dofs; ++i)
          {
            if(dof == glob_dofs[i])
            {
              loc1_dofs[i] = loc;
              return i;
            }
          }
          ASSERT(num_dofs < max_);
          glob_dofs[num_dofs] = dof;
          loc1_dofs[num_dofs] = loc;
          loc2_dofs[num_dofs] = -1;
          return num_dofs++;
        }

        int push_2(Index dof, int loc)
        {
          for(int i(0); i < num_dofs; ++i)
          {
            if(dof == glob_dofs[i])
            {
              loc2_dofs[i] = loc;
              return i;
            }
          }
          ASSERT(num_dofs < max_);
          glob_dofs[num_dofs] = dof;
          loc1_dofs[num_dofs] = -1;
          loc2_dofs[num_dofs] = loc;
          return num_dofs++;
        }

        int loc_1(int k) const
        {
          return loc1_dofs[k];
        }

        int loc_2(int k) const
        {
          return loc2_dofs[k];
        }

        int get_num_local_dofs() const
        {
          return num_dofs;
        }

        Index get_index(int i) const
        {
          return glob_dofs[i];
        }
      }; // class CommonDofMap<...>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Basic assembly task base class for a single finite element space with pairwise assembly support
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Space_, TrafoTags trafo_config_, TrafoTags facet_trafo_config_, SpaceTags space_config_>
    class TraceAssemblyJumpTaskBase1
    {
    public:
      // this task supports pairwise assembly
      static constexpr bool assemble_pairwise = true;
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
      TrafoEvaluator trafo_eval_1, trafo_eval_2;
      /// the trafo facet evaluator
      TrafoFacetEvaluator trafo_facet_eval;
      /// the space evaluator
      SpaceEvaluator space_eval_1, space_eval_2;
      /// the trafo evaluation data
      TrafoEvalData trafo_data_1, trafo_data_2;
      /// the trafo facet evaluation data
      TrafoFacetEvalData trafo_facet_data;
      /// the space evaluation data
      SpaceEvalData space_data_1, space_data_2;
      /// null-basis data
      typename AsmTraits::BasisData null_basis_data;
      /// the space dof-mapping
      typename AsmTraits::DofMapping dof_mapping_1, dof_mapping_2;
      // common DOF mapping
      static constexpr int max_common_dofs = 2 * AsmTraits::max_local_test_dofs;
      Intern::CommonDofMap<max_common_dofs> common_map;

      /// are assembling a cell pair
      bool cell_pair;

      /// the internal cell facet orientation code
      int cell_facet_ori_1, cell_facet_ori_2;

      /// local facet trafo matrices and vectors
      Tiny::Matrix<DataType, shape_dim, facet_dim> face_mat_1, face_mat_2;
      Tiny::Matrix<DataType, facet_dim, facet_dim> ori_mat_1, ori_mat_2;
      Tiny::Vector<DataType, shape_dim> face_vec_1, face_vec_2;
      Tiny::Vector<DataType, facet_dim> ori_vec_1, ori_vec_2;

      /// current cubature point on reference cell
      Tiny::Vector<DataType, shape_dim> cur_point_1, cur_point_2;
      /// current cubature point on reference facet
      Tiny::Vector<DataType, facet_dim> cur_point_facet;

    public:
      explicit TraceAssemblyJumpTaskBase1(const SpaceType& space_) :
        space(space_),
        trafo(space.get_trafo()),
        trafo_eval_1(trafo),
        trafo_eval_2(trafo),
        trafo_facet_eval(trafo),
        space_eval_1(space),
        space_eval_2(space),
        trafo_data_1(),
        trafo_data_2(),
        trafo_facet_data(),
        space_data_1(),
        space_data_2(),
        dof_mapping_1(space),
        dof_mapping_2(space),
        cell_pair(false),
        cell_facet_ori_1(0),
        cell_facet_ori_2(0)
      {
        null_basis_data.format();
        face_mat_1.format();
        face_mat_2.format();
        ori_mat_1.format();
        ori_mat_2.format();
        face_vec_1.format();
        face_vec_2.format();
        ori_vec_1.format();
        ori_vec_2.format();
      }

      virtual ~TraceAssemblyJumpTaskBase1() = default;

      void prepare(Index facet, Index cell, int local_facet, int facet_ori)
      {
        // we prepare for a single cell here
        cell_pair = false;

        // compute facet trafo
        Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_1, face_vec_1, local_facet);

        // compute orientation trafo
        Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_1, ori_vec_1, facet_ori);

        // compute orientation of actual cell facet
        cell_facet_ori_1 = Geometry::Intern::CongruencySampler<FacetType>::orientation(facet_ori)
          * Shape::ReferenceCell<ShapeType>::facet_orientation(local_facet);

        // prepare evaluators
        trafo_facet_eval.prepare(facet);
        trafo_eval_1.prepare(cell);
        space_eval_1.prepare(trafo_eval_1);

        // prepare dof mapping
        dof_mapping_1.prepare(cell);

        // build local dof map
        common_map.clear();
        for(int i(0); i < dof_mapping_1.get_num_local_dofs(); ++i)
          common_map.push_1(dof_mapping_1.get_index(i), i);
      }

      void prepare(Index facet, Index cell_1, Index cell_2, int local_facet_1, int local_facet_2, int facet_ori_1, int facet_ori_2)
      {
        // we prepare for a cell pair here
        cell_pair = true;

        // compute facet trafo
        Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_1, face_vec_1, local_facet_1);
        Geometry::Intern::FaceRefTrafo<ShapeType, facet_dim>::compute(face_mat_2, face_vec_2, local_facet_2);

        // compute orientation trafo
        Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_1, ori_vec_1, facet_ori_1);
        Geometry::Intern::CongruencyTrafo<FacetType>::compute(ori_mat_2, ori_vec_2, facet_ori_2);

        // compute orientation of actual cell facet
        cell_facet_ori_1 = Geometry::Intern::CongruencySampler<FacetType>::orientation(facet_ori_1)
          * Shape::ReferenceCell<ShapeType>::facet_orientation(local_facet_1);
        cell_facet_ori_2 = Geometry::Intern::CongruencySampler<FacetType>::orientation(facet_ori_2)
          * Shape::ReferenceCell<ShapeType>::facet_orientation(local_facet_2);

        // prepare evaluators
        trafo_facet_eval.prepare(facet);
        trafo_eval_1.prepare(cell_1);
        trafo_eval_2.prepare(cell_2);
        space_eval_1.prepare(trafo_eval_1);
        space_eval_2.prepare(trafo_eval_2);

        // prepare dof mapping
        dof_mapping_1.prepare(cell_1);
        dof_mapping_2.prepare(cell_2);

        // build local dof map
        common_map.clear();
        for(int i(0); i < dof_mapping_1.get_num_local_dofs(); ++i)
          common_map.push_1(dof_mapping_1.get_index(i), i);
        for(int i(0); i < dof_mapping_2.get_num_local_dofs(); ++i)
          common_map.push_2(dof_mapping_2.get_index(i), i);
      }

      void prepare_point(Tiny::Vector<DataType, facet_dim>& pt)
      {
        // set cubature point
        cur_point_facet = pt;

        // compute facet trafo data
        trafo_facet_eval(trafo_facet_data, cur_point_facet);

        // compute normal vector on facet trafo data
        trafo_facet_data.normal = Tiny::orthogonal(trafo_facet_data.jac_mat).normalize();

        // transform point to local facet on reference cell
        cur_point_1 = (face_mat_1 * ((ori_mat_1 * cur_point_facet) + ori_vec_1)) + face_vec_1;

        // compute trafo data
        trafo_eval_1(trafo_data_1, cur_point_1);

        // compute normal vector
        trafo_data_1.normal = trafo_facet_data.normal;

        // ensure that the normal vector "points outside"
        if(cell_facet_ori_1 < 0)
          trafo_data_1.normal.negate();

        // compute space data
        space_eval_1(space_data_1, trafo_data_1);

        if(!cell_pair)
          return;

        // transform point to local facet on reference cell
        cur_point_2 = (face_mat_2 * ((ori_mat_2 * cur_point_facet) + ori_vec_2)) + face_vec_2;

        // compute trafo data
        trafo_eval_2(trafo_data_2, cur_point_2);

        // compute normal vector
        trafo_data_2.normal = trafo_facet_data.normal;

        // ensure that the normal vector "points outside"
        if(cell_facet_ori_2 < 0)
          trafo_data_2.normal.negate();

        // compute space data
        space_eval_2(space_data_2, trafo_data_2);
      }

      void finish()
      {
        if(cell_pair)
        {
          dof_mapping_2.finish();
          space_eval_2.finish();
          trafo_eval_2.finish();
        }
        dof_mapping_1.finish();
        space_eval_1.finish();
        trafo_eval_1.finish();
        trafo_facet_eval.finish();
      }
    }; // class TraceAssemblyJumpTaskBase1<...>

    /**
     * \brief Jump Matrix assembly task CRTP base-class for identical test-/trial-space
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
    class TraceAssemblyJumpMatrixTaskCRTP1 :
      public TraceAssemblyJumpTaskBase1<typename Matrix_::DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_>
    {
    public:
      /// our base class
      typedef TraceAssemblyJumpTaskBase1<typename Matrix_::DataType, Space_, trafo_config_, facet_trafo_config_ | TrafoTags::jac_det, space_config_> BaseClass;

      /// the data-type of the matrix
      typedef typename Matrix_::DataType DataType;
      /// the value-type of the matrix
      typedef typename Matrix_::ValueType ValueType;

      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task has no combine
      static constexpr bool need_combine = false;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      using typename BaseClass::SpaceBasisData;

      using BaseClass::max_common_dofs;

    protected:
      /// the matrix that is to be assembled
      Matrix_& matrix;
      /// the cubature rule used for integration
      typename BaseClass::CubatureRuleType cubature_rule;
      /// the local matrix to be assembled
      Tiny::Matrix<ValueType, max_common_dofs, max_common_dofs> local_matrix;
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
      explicit TraceAssemblyJumpMatrixTaskCRTP1(Matrix_& matrix_, const Space_& space_,
        const Cubature::DynamicFactory& cubature_factory, DataType alpha_) :
        BaseClass(space_),
        matrix(matrix_),
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
       * \param[in] tau_f
       * The \transient facet trafo evaluation data for the current cubature point.
       *
       * \param[in] is_pair
       * Specifies whether the point is set for a pair of cells
       *
       * \param[in] tau_1
       * The \transient trafo evaluation data for the current cubature point on the first cell.
       *
       * \param[in] tau_2
       * The \transient trafo evaluation data for the current cubature point on the second cell.
       */
      void set_point(TrafoFacetEvalData& tau_f, bool is_pair, TrafoEvalData& tau_1, TrafoEvalData& tau_2);

      /**
       * \brief Evaluates the assembly for the a test-/trial-basis function pair
       *
       * \param[inout] val
       * A \transient reference to the value that is to be updated
       *
       * \param[in] weight
       * The current cubature weight
       *
       * \param[in] is_pair
       * Specifies whether the point is set for a pair of cells
       *
       * \param[in] psi_1
       * The \transient evaluation data of the current test basis function on the first cell.
       *
       * \param[in] psi_2
       * The \transient evaluation data of the current test basis function on the second cell.
       *
       * \param[in] phi_1
       * The \transient evaluation data of the current trial basis function on the first cell.
       *
       * \param[in] phi_2
       * The \transient evaluation data of the current trial basis function on the second cell.
       */
      void eval(ValueType& val, const DataType& weight, bool is_pair,
        const SpaceBasisData& psi_1, const SpaceBasisData& psi_2,
        const SpaceBasisData& phi_1, const SpaceBasisData& phi_2);
#endif // DOXYGEN

      /**
       * \brief Performs the local assembly.
       */
      void assemble()
      {
        // format local matrix
        local_matrix.format();

        // fetch number of local dofs
        const int num_loc_dofs = this->common_map.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // prepare point
          static_cast<Derived_&>(*this).prepare_point(cubature_rule.get_point(k));

          // evaluate operator
          static_cast<Derived_&>(*this).set_point(this->trafo_facet_data, this->cell_pair, this->trafo_data_1, this->trafo_data_2);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            const int i_1 = this->common_map.loc_1(i);
            const int i_2 = this->common_map.loc_2(i);
            const SpaceBasisData& psi_1 = (i_1 > -1 ? this->space_data_1.phi[i_1] : this->null_basis_data);
            const SpaceBasisData& psi_2 = (i_2 > -1 ? this->space_data_2.phi[i_2] : this->null_basis_data);

            // trial function loop
            for(int j(0); j < num_loc_dofs; ++j)
            {
              const int j_1 = this->common_map.loc_1(j);
              const int j_2 = this->common_map.loc_2(j);
              const SpaceBasisData& phi_1 = (j_1 > -1 ? this->space_data_1.phi[j_1] : this->null_basis_data);
              const SpaceBasisData& phi_2 = (j_2 > -1 ? this->space_data_2.phi[j_2] : this->null_basis_data);

              // evaluate operator and integrate
              static_cast<Derived_&>(*this).eval(local_matrix(i,j),
                this->trafo_facet_data.jac_det * cubature_rule.get_weight(k),
                this->cell_pair, psi_1, psi_2, phi_1, phi_2);
            } // continue with next trial function
          } // continue with next test function
        } // continue with next cubature point
      }

      /**
       * \brief Scatters the local assembly.
       */
      void scatter()
      {
        scatter_axpy(local_matrix, this->common_map, this->common_map, scatter_alpha);
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class TraceAssemblyJumpMatrixTaskCRTP1<...>

    /**
     * \brief Matrix assembly job for the edge-oriented jump stabilization operator
     *
     * This class implements the TraceAssemblyJob interface to assemble the jump stabilization operator:
     *   \f[J(\varphi,\psi) = \gamma \sum_E (s\cdot J_E)^{p} \int_E [\nabla \varphi]\cdot[\nabla\psi]\f]
     *
     * \tparam Matrix_
     * The matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Space_>
    class TraceAssemblyJumpStabilizationMatrixJob
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      class Task :
        public TraceAssemblyJumpMatrixTaskCRTP1<Task, Matrix_, Space_, TrafoTags::none, TrafoTags::jac_det, SpaceTags::grad>
      {
      protected:
        /// our base-class typedef
        typedef TraceAssemblyJumpMatrixTaskCRTP1<Task, Matrix_, Space_, TrafoTags::none, TrafoTags::jac_det, SpaceTags::grad> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        using typename BaseClass::TrafoFacetEvalData;
        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::SpaceBasisData;

        /// stabilization parameters
        const DataType jacdet_scale, jacdet_power;

        /// h-dependent point-wise scaling factor
        DataType point_scale;

      public:
        explicit Task(TraceAssemblyJumpStabilizationMatrixJob& job) :
          BaseClass(job.matrix, job.space, job.cubature_factory, job.gamma),
          jacdet_scale(job.jacdet_scale),
          jacdet_power(job.jacdet_power),
          point_scale()
        {
        }

        void set_point(TrafoFacetEvalData& tau_f, bool DOXY(is_pair), TrafoEvalData& DOXY(tau_1), TrafoEvalData& DOXY(tau_2))
        {
          // pre-compute scaling factor
          point_scale = Math::pow(jacdet_scale * tau_f.jac_det, jacdet_power);
        }

        void eval(ValueType& val, const DataType& weight, bool DOXY(is_pair),
          const SpaceBasisData& psi_1, const SpaceBasisData& psi_2,
          const SpaceBasisData& phi_1, const SpaceBasisData& phi_2)
        {
          // integrate < [psi], [phi] >
          Tiny::add_id(val, weight * point_scale * Tiny::dot(psi_2.grad - psi_1.grad, phi_2.grad - phi_1.grad));
        }
      }; // class Task

    protected:
      /// a reference to the matrix that is to be assembled
      Matrix_& matrix;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;
      /// the scaling factor for the assembly.
      DataType gamma;
      /// stabilization parameters
      DataType jacdet_scale, jacdet_power;

    public:
      /**
       * \brief Constructor
       *
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] space_
       * A \resident reference to the (test) space to be used for the discretization.
       *
       * \param[in] cubature_
       * The name of the cubature rule that is to be used for integration.
       *
       * \param[in] gamma_
       * The scaling factor for the entire operator.
       *
       * \param[in] jacdet_scale_
       * The scaling factor for the Jacobian determinant factor, defaults to 2.
       *
       * \param[in] jacdet_power_
       * The exponential power for the Jacobian determinant factor, defaults to 2.
       */
      explicit TraceAssemblyJumpStabilizationMatrixJob(Matrix_& matrix_, const Space_& space_, String cubature_,
        DataType gamma_ = DataType(1), DataType jacdet_scale_ = DataType(2), DataType jacdet_power_ = DataType(2)) :
        matrix(matrix_),
        space(space_),
        cubature_factory(cubature_),
        gamma(gamma_),
        jacdet_scale(jacdet_scale_),
        jacdet_power(jacdet_power_)
      {
      }
    }; // class TraceAssemblyJumpStabilizationMatrixJob<...>

    /**
     * \brief Assembles edge-oriented jump stabilization operator into a matrix
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
     *
     * \param[inout] matrix
     * The \transient matrix that is to be assembled.
     *
     * \param[in] space
     * A \transient reference to the finite-element to be used as the test- and trial-space.
     *
     * \param[in] cubature
     * The name of the cubature rule that is to be used for the assembly.
     *
     * \param[in] gamma
     * The scaling factor for the entire operator.
     *
     * \param[in] jacdet_scale
     * The scaling factor for the Jacobian determinant factor, defaults to 2.
     *
     * \param[in] jacdet_power
     * The exponential power for the Jacobian determinant factor, defaults to 2.
     */
    template<typename Trafo_, typename Matrix_, typename Space_>
    void assemble_jump_stabilization_matrix(TraceAssembler<Trafo_>& trace_asm,
      Matrix_& matrix, const Space_& space, const String& cubature,
      const typename Matrix_::DataType gamma = typename Matrix_::DataType(1),
      const typename Matrix_::DataType jacdet_scale = typename Matrix_::DataType(2),
      const typename Matrix_::DataType jacdet_power = typename Matrix_::DataType(2))
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "trace assembler and space have different trafos");
      XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix row count");

      TraceAssemblyJumpStabilizationMatrixJob<Matrix_, Space_> job(matrix, space, cubature,
        gamma, jacdet_scale, jacdet_power);
      trace_asm.assemble(job);
    }

    /**
     * \brief Matrix assembly job for the mass jump operator
     *
     * This class implements the TraceAssemblyJob interface to assemble the mass jump operator:
     *   \f[J(\varphi,\psi) = \alpha \sum_E \int_E [\varphi]\cdot[\psi]\f]
     *
     * \tparam Matrix_
     * The matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used as test- and trial-space.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Space_>
    class TraceAssemblyJumpMassMatrixJob
    {
    public:
      typedef typename Matrix_::DataType DataType;
      typedef typename Matrix_::ValueType ValueType;

      class Task :
        public TraceAssemblyJumpMatrixTaskCRTP1<Task, Matrix_, Space_, TrafoTags::none, TrafoTags::none, SpaceTags::value>
      {
      protected:
        /// our base-class typedef
        typedef TraceAssemblyJumpMatrixTaskCRTP1<Task, Matrix_, Space_, TrafoTags::none, TrafoTags::none, SpaceTags::value> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        using typename BaseClass::TrafoFacetEvalData;
        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::SpaceBasisData;

      public:
        explicit Task(TraceAssemblyJumpMassMatrixJob& job) :
          BaseClass(job.matrix, job.space, job.cubature_factory, job.alpha)
        {
        }

        void eval(ValueType& val, const DataType& weight, bool,
          const SpaceBasisData& psi_1, const SpaceBasisData& psi_2,
          const SpaceBasisData& phi_1, const SpaceBasisData& phi_2)
        {
          Tiny::add_id(val, weight * (psi_2.value - psi_1.value) * (phi_2.value - phi_1.value));
        }
      }; // class Task

    protected:
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
       * \param[inout] matrix_
       * A \resident reference to the matrix that is to be assembled.
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
      explicit TraceAssemblyJumpMassMatrixJob(Matrix_& matrix_, const Space_& space_, String cubature_, DataType alpha_ = DataType(1)) :
        matrix(matrix_),
        space(space_),
        cubature_factory(cubature_),
        alpha(alpha_)
      {
      }
    }; // class TraceAssemblyJumpMassMatrixJob<...>

    /**
     * \brief Assembles mass jump operator into a matrix
     *
     * \param[inout] trace_asm
     * A \transient reference to the trace assembler that is to be used for the assembly.
     *
     * \param[inout] matrix
     * The \transient matrix that is to be assembled.
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
    template<typename Trafo_, typename Matrix_, typename Space_>
    void assemble_jump_mass_matrix(TraceAssembler<Trafo_>& trace_asm,
      Matrix_& matrix, const Space_& space, const String& cubature,
      const typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
    {
      XASSERTM(trace_asm.get_trafo() == space.get_trafo(), "trace assembler and space have different trafos");
      XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix column count");
      XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix row count");

      TraceAssemblyJumpMassMatrixJob<Matrix_, Space_> job(matrix, space, cubature, alpha);
      trace_asm.assemble(job);
    }
  } // namespace Assembly
} // namespace FEAT
