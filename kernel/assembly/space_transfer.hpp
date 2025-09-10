// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/geometry/mesh_permutation.hpp>



namespace FEAT
{
  namespace Assembly
  {
    class SpaceTransfer
    {
    public:
      /**
       * \brief Exception for singular mass matrix errors
       *
       * This exception is thrown if the inversion of a local mass matrix fails, thus resulting
       * in invalid entries in a prolongation/truncation/transfer matrix.
       */
      class LocalMassMatrixSingularException :
        public FEAT::Exception
      {
      public:
        LocalMassMatrixSingularException() :
          Exception("local mass matrix inversion error")
        {
        }
      }; // class LocalMassMatrixSingularException

      /**
       * \brief Assembles an inter-space transfer matrix and its corresponding weight vector.
       *
       * To obtain the final transfer matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A \transient reference to the inter-space transfer matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the transfer matrix.
       *
       * \param[in] target_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] source_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_>
      static void assemble_transfer(
        Matrix_& matrix,
        Vector_& vector,
        const TargetSpace_& target_space,
        const SourceSpace_& source_space,
        const String& cubature_name)
      {
        // source and target space must be defined on the same trafo and therefore on the same mesh
        XASSERTM(&target_space.get_trafo() == &source_space.get_trafo(), "source and target space must be defined on the same trafo");

        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == target_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == source_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector.size() == target_space.get_num_dofs(), "invalid vector size");

        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename TargetSpace_::TrafoType TrafoType;
        typedef typename SourceSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename TargetSpace_::DofMappingType TargetDofMapping;
        typedef typename SourceSpace_::DofMappingType SourceDofMapping;

        // typedefs for trafo evaluators
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

        // typedefs for space evaluators
        typedef typename TargetSpace_::template Evaluator<TrafoEvaluator>::Type TargetSpaceEvaluator;
        typedef typename SourceSpace_::template Evaluator<TrafoEvaluator>::Type SourceSpaceEvaluator;

        // define fine and coarse mesh trafo configurations
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value> TargetSpaceConfigTraits;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value> SourceSpaceConfigTraits;
        static constexpr TrafoTags trafo_config = TrafoTags::jac_det | TargetSpaceConfigTraits::trafo_config | SourceSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;
        TrafoEvalData trafo_data;

        // typedef for space data
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType TargetSpaceEvalData;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType SourceSpaceEvalData;
        TargetSpaceEvalData target_space_data;
        SourceSpaceEvalData source_space_data;

        // create matrix scatter-axpy
        typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
        typename Vector_::ScatterAxpy scatter_vaxpy(vector);

        // create DOF-mappings
        TargetDofMapping target_dof_mapping(target_space);
        SourceDofMapping source_dof_mapping(source_space);

        // fetch the trafos
        const TrafoType& trafo = target_space.get_trafo();

        // create trafo evaluators
        TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        TargetSpaceEvaluator target_space_eval(target_space);
        SourceSpaceEvaluator source_space_eval(source_space);

        // create the cubature rules
        Cubature::DynamicFactory cubature_factory(cubature_name);
        typename Intern::CubatureTraits<TrafoEvaluator>::RuleType cubature(Cubature::ctor_factory, cubature_factory);

        // allocate fine-mesh mass matrix
        Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, TargetSpaceEvaluator::max_local_dofs> mass;

        // allocate local matrix data for inter-space mass matrix
        Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, SourceSpaceEvaluator::max_local_dofs> lmd, lid;

        // allocate local vector data for weight vector
        Tiny::Vector<DataType, TargetSpaceEvaluator::max_local_dofs> lvd;

        // pivot array for factorization
        int pivot[TargetSpaceEvaluator::max_local_dofs];

        // loop over all coarse mesh cells
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          source_space_eval.prepare(trafo_eval);
          target_space_eval.prepare(trafo_eval);

          // fetch number of local DOFs
          const int source_num_loc_dofs = source_space_eval.get_num_local_dofs();
          const int target_num_loc_dofs = target_space_eval.get_num_local_dofs();

          // prepare dof-mapping
          source_dof_mapping.prepare(cell);
          target_dof_mapping.prepare(cell);

          // format local matrices
          mass.format();
          lmd.format();
          lvd.format();

          // loop over all cubature points and integrate
          for(int k(0); k < cubature.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature.get_point(k));

            // compute basis function data
            target_space_eval(target_space_data, trafo_data);
            source_space_eval(source_space_data, trafo_data);

            // pre-compute cubature weight
            const DataType omega = trafo_data.jac_det * cubature.get_weight(k);

            // target mesh test function loop
            for(int i(0); i < target_num_loc_dofs; ++i)
            {
              // target mesh trial function loop
              for(int j(0); j < target_num_loc_dofs; ++j)
              {
                mass(i,j) += omega * target_space_data.phi[i].value * target_space_data.phi[j].value;
              }

              // source mesh trial function loop
              for(int j(0); j < source_num_loc_dofs; ++j)
              {
                lmd(i,j) += omega * target_space_data.phi[i].value * source_space_data.phi[j].value;
              }
            }
          }

          // invert target space mass matrix
          Math::invert_matrix(target_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

          // Note:
          // Usually, one would check whether the determinant returned by the invert_matrix
          // function is normal. However, this can lead to false alerts when assembling in
          // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
          // the determinant can become subnormal or even (numerically) zero although the
          // condition number of the matrix is still fine and the inversion was successful.
          // Therefore, we first multiply the (hopefully) inverted mass matrix by the
          // inter-level mass matrix and check whether the Frobenius norm of the result
          // is normal. If our matrix inversion failed, the result is virtually guaranteed
          // to be garbage, so this should serve well enough as a sanity check.

          // compute X := M^{-1}*N
          lid.set_mat_mat_mult(mass, lmd);

          // sanity check for matrix inversion
          if(!Math::isnormal(lid.norm_frobenius()))
            throw LocalMassMatrixSingularException();

          // incorporate local matrix
          scatter_maxpy(lid, target_dof_mapping, source_dof_mapping);

          // update weights
          lvd.format(DataType(1));
          scatter_vaxpy(lvd, target_dof_mapping);

          // finish dof-mappings
          source_dof_mapping.finish();
          target_dof_mapping.finish();

          // finish evaluators
          target_space_eval.finish();
          source_space_eval.finish();
          trafo_eval.finish();

          // go for next coarse mesh cell
        }
      }

      /**
       * \brief Assembles an inter-space transfer matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble prolongation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A \transient reference to the transfer matrix that is to be assembled.
       *
       * \param[in] target_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] source_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename TargetSpace_,
        typename SourceSpace_>
      static void assemble_transfer_direct(
        Matrix_& matrix,
        const TargetSpace_& target_space,
        const SourceSpace_& source_space,
        const String& cubature_name)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        assemble_transfer(matrix, weight, target_space, source_space, cubature_name);

        // scale prolongation matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);
      }

      /**
       * \brief Transfers a primal vector and assembles a compatible weight vector
       *
       * To obtain the final transferred vector, one needs to invert the weight vector component-wise
       * and scale the transferred vector component-wise by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c component_product
       * operations of the vector container, resp.
       *
       * \attention
       * In the case of global vectors, the weight vector has to be synchronized via sync_0 before
       * the component-wise inversion and scaling to obtain the correctly scaled primal vector.
       *
       * \param[in,out] vector_f
       * A \transient reference to the fine-mesh vector that is to be assembled
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in,out] vector_w
       * A \transient reference to the fine-mesh weight vector that is to be assembled.
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_c
       * A \transient reference to the coarse-mesh vector that is to be transferred.
       *
       * \param[in] target_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] source_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_>
      static void transfer_vector(
        Vector_& vector_f,
        Vector_& vector_w,
        const Vector_& vector_c,
        const TargetSpace_& target_space,
        const SourceSpace_& source_space,
        const String& cubature_name)
      {
        // source and target space must be defined on the same trafo and therefore on the same mesh
        XASSERTM(&target_space.get_trafo() == &source_space.get_trafo(), "source and target space must be defined on the same trafo");

        // validate matrix and vector dimensions
        XASSERTM(vector_f.size() == vector_w.size(), "invalid vector sizes");
        XASSERTM(vector_f.size() == target_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector_c.size() == source_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector_f.size() == target_space.get_num_dofs(), "invalid vector size");

        typedef typename Vector_::DataType DataType;
        typedef typename Vector_::ValueType ValueType;

        // typedefs for trafos, mesh and shape
        typedef typename TargetSpace_::TrafoType TrafoType;
        typedef typename SourceSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename TargetSpace_::DofMappingType TargetDofMapping;
        typedef typename SourceSpace_::DofMappingType SourceDofMapping;

        // typedefs for trafo evaluators
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

        // typedefs for space evaluators
        typedef typename TargetSpace_::template Evaluator<TrafoEvaluator>::Type TargetSpaceEvaluator;
        typedef typename SourceSpace_::template Evaluator<TrafoEvaluator>::Type SourceSpaceEvaluator;

        // define target and source mesh trafo configurations
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value> TargetSpaceConfigTraits;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value> SourceSpaceConfigTraits;
        static constexpr TrafoTags target_trafo_config = TrafoTags::jac_det | TargetSpaceConfigTraits::trafo_config | SourceSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename TrafoEvaluator::template ConfigTraits<target_trafo_config>::EvalDataType TrafoEvalData;
        TrafoEvalData trafo_data;

        // typedef for space data
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType TargetSpaceEvalData;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType SourceSpaceEvalData;
        TargetSpaceEvalData target_space_data;
        SourceSpaceEvalData source_space_data;

        // typedefs for cubature rule and factory
        typedef typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType CubatureRuleType;

        // create gather/scatter-axpys
        typename Vector_::GatherAxpy gather_c(vector_c);
        typename Vector_::ScatterAxpy scatter_f(vector_f);
        typename Vector_::ScatterAxpy scatter_w(vector_w);

        // create DOF-mappings
        TargetDofMapping target_dof_mapping(target_space);
        SourceDofMapping source_dof_mapping(source_space);

        // fetch the trafos
        const TrafoType& trafo = target_space.get_trafo();

        // create trafo evaluators
        TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        TargetSpaceEvaluator target_space_eval(target_space);
        SourceSpaceEvaluator source_space_eval(source_space);

        // create the cubature rules
        Cubature::DynamicFactory cubature_factory(cubature_name);
        CubatureRuleType cubature(Cubature::ctor_factory, cubature_factory);

        // allocate fine-mesh mass matrix
        Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, TargetSpaceEvaluator::max_local_dofs> mass;

        // allocate local matrix data for interlevel-mesh mass matrix
        Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, SourceSpaceEvaluator::max_local_dofs> lmd, lid;

        // allocate local vector data for weight/coarse/fine vector
        Tiny::Vector<ValueType, TargetSpaceEvaluator::max_local_dofs, TargetSpaceEvaluator::max_local_dofs> lv_w, lv_f;
        Tiny::Vector<ValueType, SourceSpaceEvaluator::max_local_dofs, SourceSpaceEvaluator::max_local_dofs> lv_c;

        // pivot array for factorization
        int pivot[TargetSpaceEvaluator::max_local_dofs];

        // loop over all coarse mesh cells
        for(Index ccell(0); ccell < trafo_eval.get_num_cells(); ++ccell)
        {
          // prepare coarse trafo evaluator
          trafo_eval.prepare(ccell);

          // prepare coarse space evaluator
          source_space_eval.prepare(trafo_eval);
          target_space_eval.prepare(trafo_eval);

          // prepare dof-mapping
          source_dof_mapping.prepare(ccell);
          target_dof_mapping.prepare(ccell);

          // prepare coarse mesh dof-mapping
          source_dof_mapping.prepare(ccell);

          // fetch number of local coarse DOFs
          const int source_num_loc_dofs = source_space_eval.get_num_local_dofs();
          const int target_num_loc_dofs = target_space_eval.get_num_local_dofs();

          // gather local coarse vector
          lv_c.format();
          gather_c(lv_c, source_dof_mapping);

          // format local matrices
          mass.format();
          lmd.format();
          lv_w.format();
          lv_f.format();

          // loop over all cubature points and integrate
          for(int k(0); k < cubature.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature.get_point(k));

            // compute basis function data
            target_space_eval(target_space_data, trafo_data);
            source_space_eval(source_space_data, trafo_data);

            // pre-compute cubature weight
            const DataType omega = trafo_data.jac_det * cubature.get_weight(k);

            // target mesh test function loop
            for(int i(0); i < target_num_loc_dofs; ++i)
            {
              // target mesh trial function loop
              for(int j(0); j < target_num_loc_dofs; ++j)
              {
                mass(i,j) += omega * target_space_data.phi[i].value * target_space_data.phi[j].value;
              }

              // source mesh trial function loop
              for(int j(0); j < source_num_loc_dofs; ++j)
              {
                lmd(i,j) += omega * target_space_data.phi[i].value * source_space_data.phi[j].value;
              }
            }
            // go for next cubature point
          }
          // finish coarse mesh evaluators
          target_space_eval.finish();
          trafo_eval.finish();

          // invert target mesh mass matrix
          Math::invert_matrix(target_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

          // Note:
          // Usually, one would check whether the determinant returned by the invert_matrix
          // function is normal. However, this can lead to false alerts when assembling in
          // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
          // the determinant can become subnormal or even (numerically) zero although the
          // condition number of the matrix is still fine and the inversion was successful.
          // Therefore, we first multiply the (hopefully) inverted mass matrix by the
          // inter-level mass matrix and check whether the Frobenius norm of the result
          // is normal. If our matrix inversion failed, the result is virtually guaranteed
          // to be garbage, so this should serve well enough as a sanity check.
          // compute X := M^{-1}*N

          lid.set_mat_mat_mult(mass, lmd);

          // sanity check for matrix inversion
          if(!Math::isnormal(lid.norm_frobenius()))
            throw LocalMassMatrixSingularException();

          // compute local target vector
          lv_f.format();
          for(int i(0); i < target_num_loc_dofs; ++i)
            for(int j(0); j < source_num_loc_dofs; ++j)
              lv_f[i] += lid[i][j] * lv_c[j];

           // prepare target mesh dof-mapping
          target_dof_mapping.prepare(ccell);

          // scatter local target vector
          scatter_f(lv_f, target_dof_mapping);

          // update weights
          lv_w.format(DataType(1));
          scatter_w(lv_w, target_dof_mapping);

          // finish target mesh dof-mapping
          target_dof_mapping.finish();

          // finish coarse mesh evaluators
          source_space_eval.finish();
          trafo_eval.finish();

          // finish coarse mesh dof-mapping
          source_dof_mapping.finish();

          // go for next coarse mesh cell
        }
      }

      /**
       * \brief Transfers a primal vector directly
       *
       * \attention
       * This function <b>must not</b> be used to transfer vectors for parallel (i.e. global)
       * simulations, as it will be scaled incorrectly due to missing weight vector synchronization!
       *
       * \param[in,out] vector_f
       * A \transient reference to the fine-mesh vector that is to be assembled
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_c
       * A \transient reference to the coarse-mesh vector that is to be transferred.
       *
       * \param[in] target_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] source_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_>
      static void transfer_vector_direct(
        Vector_& vector_f,
        const Vector_& vector_c,
        const TargetSpace_& target_space,
        const SourceSpace_& source_space,
        const String& cubature_name)
      {
        Vector_ vector_w = vector_f.clone(LAFEM::CloneMode::Layout);
        vector_w.format();

        transfer_vector(vector_f, vector_w, vector_c, target_space, source_space, cubature_name);

        // finally, scale target mesh vector by inverse weights
        vector_w.component_invert(vector_w);
        vector_f.component_product(vector_f, vector_w);
      }
    }; // class SpaceTransfer
  } // namespace Assembly
} // namespace FEAT
