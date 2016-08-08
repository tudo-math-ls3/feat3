#pragma once
#ifndef KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
#define KERNEL_ASSEMBLY_GRID_TRANSFER_HPP 1

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Grid-Transfer assembly class template
     *
     * This class template implements the assembly of grid transfer operators.
     *
     * \author Peter Zajac
     */
    class GridTransfer
    {
    public:
      /**
       * \brief Assembles a prolongation matrix and its corresponding weight vector.
       *
       * To obtain the final prolongation matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A reference to the prolongation matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A reference to the weight vector for the prolongation matrix.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] fine_space
       * A reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A reference to the coarse-mesh trial-space to be used.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_,
        typename CubatureFactory_>
      static void assemble_prolongation(
        Matrix_& matrix,
        Vector_& vector,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const CubatureFactory_& cubature_factory)
      {
        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename FineSpace_::TrafoType FineTrafoType;
        typedef typename CoarseSpace_::TrafoType CoarseTrafoType;
        typedef typename CoarseSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename FineSpace_::DofMappingType FineDofMapping;
        typedef typename CoarseSpace_::DofMappingType CoarseDofMapping;

        // typedefs for trafo evaluators
        typedef typename FineTrafoType::template Evaluator<ShapeType, DataType>::Type FineTrafoEvaluator;
        typedef typename CoarseTrafoType::template Evaluator<ShapeType, DataType>::Type CoarseTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename FineSpace_::template Evaluator<FineTrafoEvaluator>::Type FineSpaceEvaluator;
        typedef typename CoarseSpace_::template Evaluator<CoarseTrafoEvaluator>::Type CoarseSpaceEvaluator;

        // define fine and coarse mesh trafo configurations
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value> FineSpaceConfigTraits;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value> CoarseSpaceConfigTraits;
        static constexpr TrafoTags fine_trafo_config = TrafoTags::jac_det | FineSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags coarse_trafo_config = TrafoTags::jac_det | CoarseSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename FineTrafoEvaluator::template ConfigTraits<fine_trafo_config>::EvalDataType FineTrafoEvalData;
        typedef typename CoarseTrafoEvaluator::template ConfigTraits<coarse_trafo_config>::EvalDataType CoarseTrafoEvalData;
        FineTrafoEvalData fine_trafo_data;
        CoarseTrafoEvalData coarse_trafo_data;

        // typedef for space data
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType FineSpaceEvalData;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType CoarseSpaceEvalData;
        FineSpaceEvalData fine_space_data;
        CoarseSpaceEvalData coarse_space_data;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;

        // create matrix scatter-axpy
        typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
        typename Vector_::ScatterAxpy scatter_vaxpy(vector);

        // create DOF-mappings
        FineDofMapping fine_dof_mapping(fine_space);
        CoarseDofMapping coarse_dof_mapping(coarse_space);

        // fetch the trafos
        const FineTrafoType& fine_trafo = fine_space.get_trafo();
        const CoarseTrafoType& coarse_trafo = coarse_space.get_trafo();

        // create trafo evaluators
        FineTrafoEvaluator fine_trafo_eval(fine_trafo);
        CoarseTrafoEvaluator coarse_trafo_eval(coarse_trafo);

        // create space evaluators
        FineSpaceEvaluator fine_space_eval(fine_space);
        CoarseSpaceEvaluator coarse_space_eval(coarse_space);

        // create the cubature rules
        CubatureRuleType coarse_cubature, fine_cubature(Cubature::ctor_factory, cubature_factory);
        Cubature::RefineFactoryCore::create(coarse_cubature, fine_cubature);

        // allocate fine-mesh mass matrix
        Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, FineSpaceEvaluator::max_local_dofs> mass;

        // allocate local matrix data for interlevel-mesh mass matrix
        Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, CoarseSpaceEvaluator::max_local_dofs> lmd, lid;

        // allocate local vector data for weight vector
        Tiny::Vector<DataType, FineSpaceEvaluator::max_local_dofs> lvd;

        // pivaot array for factorisation
        int pivot[FineSpaceEvaluator::max_local_dofs];

        // calculate child count
        const Index num_children = fine_trafo_eval.get_num_cells() / coarse_trafo_eval.get_num_cells();

        // loop over all coarse mesh cells
        for(Index ccell(0); ccell < coarse_trafo_eval.get_num_cells(); ++ccell)
        {
          // prepare coarse trafo evaluator
          coarse_trafo_eval.prepare(ccell);

          // prepare coarse space evaluator
          coarse_space_eval.prepare(coarse_trafo_eval);

          // fetch number of local coarse DOFs
          int coarse_num_loc_dofs = coarse_space_eval.get_num_local_dofs();

          // prepare coarse mesh dof-mapping
          coarse_dof_mapping.prepare(ccell);

          // loop over all child cells
          for(Index child(0); child < num_children; ++child)
          {
            // calculate fine mesh cell index
            Index fcell = ccell*num_children + child;

            // prepare fine trafo evaluator
            fine_trafo_eval.prepare(fcell);

            // prepare fine space evaluator
            fine_space_eval.prepare(fine_trafo_eval);

            // fetch number of local fine DOFs
            int fine_num_loc_dofs = fine_space_eval.get_num_local_dofs();

            // format local matrices
            mass.format();
            lmd.format();
            lvd.format();

            // loop over all cubature points and integrate
            for(int k(0); k < fine_cubature.get_num_points(); ++k)
            {
              // compute coarse mesh cubature point index
              int l(int(child) * fine_cubature.get_num_points() + k);

              // compute trafo data
              fine_trafo_eval(fine_trafo_data, fine_cubature.get_point(k));
              coarse_trafo_eval(coarse_trafo_data, coarse_cubature.get_point(l));

              // compute basis function data
              fine_space_eval(fine_space_data, fine_trafo_data);
              coarse_space_eval(coarse_space_data, coarse_trafo_data);

              // fine mesh test function loop
              for(int i(0); i < fine_num_loc_dofs; ++i)
              {
                // fine mesh trial function loop
                for(int j(0); j < fine_num_loc_dofs; ++j)
                {
                  mass(i,j) += fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.phi[i].value * fine_space_data.phi[j].value;
                  // go for next fine mesh trial DOF
                }

                // coarse mesh trial function loop
                for(int j(0); j < coarse_num_loc_dofs; ++j)
                {
                  lmd(i,j) +=
                    fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.phi[i].value * coarse_space_data.phi[j].value;
                  // go for next fine mesh trial DOF
                }
                // go for next fine mesh test DOF
              }
              // go for next cubature point
            }

            // finish coarse mesh evaluators
            fine_space_eval.finish();
            fine_trafo_eval.finish();

            // invert fine mesh mass matrix
            Math::invert_matrix(fine_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

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
            {
              throw InternalError(__func__,__FILE__,__LINE__,"Local Mass Matrix inversion failed!");
            }

            // prepare fine mesh dof-mapping
            fine_dof_mapping.prepare(fcell);

            // incorporate local matrix
            scatter_maxpy(lid, fine_dof_mapping, coarse_dof_mapping);

            // update weights
            lvd.format(DataType(1));
            scatter_vaxpy(lvd, fine_dof_mapping);

            // finish fine mesh dof-mapping
            fine_dof_mapping.finish();

            // go for next child cell
          }

          // finish coarse mesh evaluators
          coarse_space_eval.finish();
          coarse_trafo_eval.finish();

          // finish coarse mesh dof-mapping
          coarse_dof_mapping.finish();

          // go for next coarse mesh cell
        }
      }

      /**
       * \brief Assembles a prolongation matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble prolongation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronisation!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A reference to the prolongation matrix that is to be assembled.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] fine_space
       * A reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A reference to the coarse-mesh trial-space to be used.
       */
      template<
        typename Matrix_,
        typename FineSpace_,
        typename CoarseSpace_,
        typename CubatureFactory_>
      static void assemble_prolongation_direct(
        Matrix_& matrix,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const CubatureFactory_& cubature_factory)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        assemble_prolongation(matrix, weight, fine_space, coarse_space, cubature_factory);

        // scale prolongation matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);
      }
    }; // class GridTransfer<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
