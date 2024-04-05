// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/assembly/asm_traits.hpp>

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
       * \param[in] space_target
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename SpaceTarget_,
        typename SpaceSource_>
      static void assemble_transfer(
        Matrix_& matrix,
        Vector_& vector,
        const SpaceTarget_& space_target,
        const SpaceSource_& space_source,
        const String& cubature_name)
      {
        // source and target space must be defined on the same trafo and therefore on the same mesh
        XASSERTM(&space_target.get_trafo() == &space_source.get_trafo(), "source and target space must be defined on the same trafo");

        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space_target.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space_source.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector.size() == space_target.get_num_dofs(), "invalid vector size");

        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename SpaceTarget_::TrafoType TrafoType;
        typedef typename SpaceSource_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename SpaceTarget_::DofMappingType TargetDofMapping;
        typedef typename SpaceSource_::DofMappingType SourceDofMapping;

        // typedefs for trafo evaluators
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;

        // typedefs for space evaluators
        typedef typename SpaceTarget_::template Evaluator<TrafoEvaluator>::Type TargetSpaceEvaluator;
        typedef typename SpaceSource_::template Evaluator<TrafoEvaluator>::Type SourceSpaceEvaluator;

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
        TargetSpaceEvalData space_target_data;
        SourceSpaceEvalData space_source_data;

        // create matrix scatter-axpy
        typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
        typename Vector_::ScatterAxpy scatter_vaxpy(vector);

        // create DOF-mappings
        TargetDofMapping target_dof_mapping(space_target);
        SourceDofMapping source_dof_mapping(space_source);

        // fetch the trafos
        const TrafoType& trafo = space_target.get_trafo();

        // create trafo evaluators
        TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        TargetSpaceEvaluator space_target_eval(space_target);
        SourceSpaceEvaluator space_source_eval(space_source);

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
          space_source_eval.prepare(trafo_eval);
          space_target_eval.prepare(trafo_eval);

          // fetch number of local DOFs
          const int source_num_loc_dofs = space_source_eval.get_num_local_dofs();
          const int target_num_loc_dofs = space_target_eval.get_num_local_dofs();

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
            space_target_eval(space_target_data, trafo_data);
            space_source_eval(space_source_data, trafo_data);

            // pre-compute cubature weight
            const DataType omega = trafo_data.jac_det * cubature.get_weight(k);

            // target mesh test function loop
            for(int i(0); i < target_num_loc_dofs; ++i)
            {
              // target mesh trial function loop
              for(int j(0); j < target_num_loc_dofs; ++j)
              {
                mass(i,j) += omega * space_target_data.phi[i].value * space_target_data.phi[j].value;
              }

              // source mesh trial function loop
              for(int j(0); j < source_num_loc_dofs; ++j)
              {
                lmd(i,j) += omega * space_target_data.phi[i].value * space_source_data.phi[j].value;
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
          space_target_eval.finish();
          space_source_eval.finish();
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
       * \param[in] space_target
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename SpaceTarget_,
        typename SpaceSource_>
      static void assemble_transfer_direct(
        Matrix_& matrix,
        const SpaceTarget_& space_target,
        const SpaceSource_& space_source,
        const String& cubature_name)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        assemble_transfer(matrix, weight, space_target, space_source, cubature_name);

        // scale prolongation matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);
      }
    }; // class SpaceTransfer
  } // namespace Assembly
} // namespace FEAT
