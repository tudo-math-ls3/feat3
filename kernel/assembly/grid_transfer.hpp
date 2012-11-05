#pragma once
#ifndef KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
#define KERNEL_ASSEMBLY_GRID_TRANSFER_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/linear_algebra.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Grid-Transfer assembly class template
     *
     * This class template implements the assembly of grid transfer operators.
     *
     * \tparam Matrix_
     * The type of the matrix that is to be assembled.
     *
     * \tparam FineSpace_
     * The type of the fine-mesh finite-element space.
     *
     * \tparam CoarseSpace_
     * The type of the coarse-mesh finite-element space.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename FineSpace_,
      typename CoarseSpace_>
    class GridTransfer
    {
    public:
      /// data type
      typedef typename Matrix_::data_type DataType;

    private:
      /// \cond internal
      struct AsmSpaceConfig :
        public Space::ConfigBase
      {
        enum
        {
          need_value = 1
        };
      };
      struct AsmTrafoConfig :
        public Trafo::ConfigBase
      {
        enum
        {
          need_jac_det = 1
        };
      };
      /// \endcond

    public:
      /**
       * \brief Assembles a prolongation matrix.
       *
       * \param[in,out] matrix
       * A reference to the prolongation matrix that is to be assembled.
       *
       * \param[in] fine_space
       * A reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       */
      static void assemble_prolongation(
        Matrix_& matrix,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {

        // typedefs for trafos, mesh and shape
        typedef typename FineSpace_::TrafoType FineTrafoType;
        typedef typename CoarseSpace_::TrafoType CoarseTrafoType;
        typedef typename CoarseTrafoType::MeshType MeshType;
        typedef typename CoarseSpace_::ShapeType ShapeType;

        // define fine and coarse mesh trafo configurations
        typedef typename CoarseSpace_::template TrafoConfig<AsmSpaceConfig> CoarseSpaceTrafoConfig;
        typedef typename FineSpace_::template TrafoConfig<AsmSpaceConfig> FineSpaceTrafoConfig;
        typedef Trafo::ConfigOr<AsmTrafoConfig, CoarseSpaceTrafoConfig> CoarseTrafoConfig;
        typedef Trafo::ConfigOr<AsmTrafoConfig, FineSpaceTrafoConfig> FineTrafoConfig;

        // typedefs for dof-mappings
        typedef typename FineSpace_::template DofMapping<>::Type FineDofMapping;
        typedef typename CoarseSpace_::template DofMapping<>::Type CoarseDofMapping;

        // typedefs for trafo evaluators and data
        typedef typename FineTrafoType::template Evaluator<ShapeType, DataType>::Type FineTrafoEvaluator;
        typedef typename CoarseTrafoType::template Evaluator<ShapeType, DataType>::Type CoarseTrafoEvaluator;
        typedef Trafo::EvalData<FineTrafoEvaluator, FineTrafoConfig> FineTrafoEvalData;
        typedef Trafo::EvalData<CoarseTrafoEvaluator, CoarseTrafoConfig> CoarseTrafoEvalData;
        FineTrafoEvalData fine_trafo_data;
        CoarseTrafoEvalData coarse_trafo_data;

        // typedefs for space evaluators and data
        typedef typename FineSpace_::template Evaluator<FineTrafoEvaluator>::Type FineSpaceEvaluator;
        typedef typename CoarseSpace_::template Evaluator<CoarseTrafoEvaluator>::Type CoarseSpaceEvaluator;
        typedef Space::EvalData<FineSpaceEvaluator, AsmSpaceConfig> FineSpaceEvalData;
        typedef Space::EvalData<CoarseSpaceEvaluator, AsmSpaceConfig> CoarseSpaceEvalData;
        FineSpaceEvalData fine_space_data;
        CoarseSpaceEvalData coarse_space_data;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;
        typedef typename Cubature::DynamicFactorySelect<CubatureRuleType>::Type CubatureFactory;

        // create matrix scatter-axpy
        MatrixScatterAxpy<Matrix_> scatter_axpy(matrix);

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
        CubatureRuleType fine_cubature(CubatureFactory::create(cubature_name));
        CubatureRuleType coarse_cubature(CubatureFactory::create("refine:" + cubature_name));

        // allocate fine-mesh mass matrix
        Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, FineSpaceEvaluator::max_local_dofs> mass;

        // allocate local matrix data for interlevel-mesh mass matrix
        LocalMatrixData<
          Tiny::Matrix<
            DataType,
            FineSpaceEvaluator::max_local_dofs,
            CoarseSpaceEvaluator::max_local_dofs>,
          FineDofMapping,
          CoarseDofMapping> lmd(fine_dof_mapping, coarse_dof_mapping);

        // pivaot array for factorisation
        Index pivot[FineSpaceEvaluator::max_local_dofs];

        // allocate weight array
        DataType* weight = new DataType[fine_space.get_num_dofs()];
        LinAlg::vec_clear(fine_space.get_num_dofs(), weight);

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
          Index coarse_num_loc_dofs = coarse_space_eval.get_num_local_dofs();

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
            Index fine_num_loc_dofs = fine_space_eval.get_num_local_dofs();

            // clear local matrices
            mass.clear();
            lmd.clear();

            // loop over all cubature points and integrate
            for(Index k(0); k < fine_cubature.get_num_points(); ++k)
            {
              // compute coarse mesh cubature point index
              Index l = child * fine_cubature.get_num_points() + k;

              // compute trafo data
              fine_trafo_data(fine_trafo_eval, fine_cubature.get_point(k));
              coarse_trafo_data(coarse_trafo_eval, coarse_cubature.get_point(l));

              // compute basis function data
              fine_space_data(fine_space_eval, fine_trafo_data);
              coarse_space_data(coarse_space_eval, coarse_trafo_data);

              // fine mesh test function loop
              for(Index i(0); i < fine_num_loc_dofs; ++i)
              {
                // fine mesh trial function loop
                for(Index j(0); j < fine_num_loc_dofs; ++j)
                {
                  mass(i,j) += fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.values[i] * fine_space_data.values[j];
                  // go for next fine mesh trial DOF
                }

                // coarse mesh trial function loop
                for(Index j(0); j < coarse_num_loc_dofs; ++j)
                {
                  lmd(i,j) +=
                    fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.values[i] * coarse_space_data.values[j];
                  // go for next fine mesh trial DOF
                }
                // go for next fine mesh test DOF
              }
              // go for next cubature point
            }

            // finish coarse mesh evaluators
            fine_space_eval.finish();
            fine_trafo_eval.finish();

            // factorise fine mesh mass matrix
            LinAlg::mat_factorise(fine_num_loc_dofs, fine_num_loc_dofs,
              Index(FineSpaceEvaluator::max_local_dofs), &mass.v[0][0], pivot);

            // solve M*X = N
            LinAlg::mat_solve_mat<false>(fine_num_loc_dofs, coarse_num_loc_dofs,
              Index(CoarseSpaceEvaluator::max_local_dofs), &lmd.v[0][0],
              Index(FineSpaceEvaluator::max_local_dofs), &mass.v[0][0], pivot);

            // prepare fine mesh dof-mapping
            fine_dof_mapping.prepare(fcell);

            // incorporate local matrix
            scatter_axpy(lmd);

            // update weights
            for(Index i(0); i < fine_dof_mapping.get_num_local_dofs(); ++i)
            {
              for(Index ic(0); ic < fine_dof_mapping.get_num_contribs(i); ++ic)
              {
                weight[fine_dof_mapping.get_index(i,ic)] += DataType(fine_dof_mapping.get_weight(i,ic));
              }
            }

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

        // finally, scale rows by weights
        for(Index i(0); i < fine_space.get_num_dofs(); ++i)
        {
          if(weight[i] > DataType(0))
            weight[i] = DataType(1) / weight[i];
        }
        Intern::RowScaler<Matrix_>::apply(matrix, weight);

        // delete weights
        delete [] weight;
      }
    }; // class GridTransfer<...>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
