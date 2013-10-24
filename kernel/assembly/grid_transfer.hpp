#pragma once
#ifndef KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
#define KERNEL_ASSEMBLY_GRID_TRANSFER_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/linear_algebra.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Matrix_>
      class RowScaler;
    }
    /// \endcond

    /**
     * \brief Grid-Transfer assembly class template
     *
     * This class template implements the assembly of grid transfer operators.
     *
     * \author Peter Zajac
     */
    class GridTransfer
    {
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
        typename CubatureFactory_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_prolongation(
        Matrix_& matrix,
        const CubatureFactory_& cubature_factory,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space)
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
        typedef typename FineSpaceEvaluator::template ConfigTraits<AsmSpaceConfig>::TrafoConfig FineSpaceTrafoConfig;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<AsmSpaceConfig>::TrafoConfig CoarseSpaceTrafoConfig;
        typedef Trafo::ConfigOr<AsmTrafoConfig, FineSpaceTrafoConfig> FineTrafoConfig;
        typedef Trafo::ConfigOr<AsmTrafoConfig, CoarseSpaceTrafoConfig> CoarseTrafoConfig;

        // typedefs for trafo data
        typedef typename FineTrafoEvaluator::template ConfigTraits<FineTrafoConfig>::EvalDataType FineTrafoEvalData;
        typedef typename CoarseTrafoEvaluator::template ConfigTraits<CoarseTrafoConfig>::EvalDataType CoarseTrafoEvalData;
        FineTrafoEvalData fine_trafo_data;
        CoarseTrafoEvalData coarse_trafo_data;

        // typedef for space data
        typedef typename FineSpaceEvaluator::template ConfigTraits<AsmSpaceConfig>::EvalDataType FineSpaceEvalData;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<AsmSpaceConfig>::EvalDataType CoarseSpaceEvalData;
        FineSpaceEvalData fine_space_data;
        CoarseSpaceEvalData coarse_space_data;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;

        // create matrix scatter-axpy
        LAFEM::ScatterAxpy<Matrix_> scatter_axpy(matrix);

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
              fine_trafo_eval(fine_trafo_data, fine_cubature.get_point(k));
              coarse_trafo_eval(coarse_trafo_data, coarse_cubature.get_point(l));

              // compute basis function data
              fine_space_eval(fine_space_data, fine_trafo_data);
              coarse_space_eval(coarse_space_data, coarse_trafo_data);

              // fine mesh test function loop
              for(Index i(0); i < fine_num_loc_dofs; ++i)
              {
                // fine mesh trial function loop
                for(Index j(0); j < fine_num_loc_dofs; ++j)
                {
                  mass(i,j) += fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.phi[i].value * fine_space_data.phi[j].value;
                  // go for next fine mesh trial DOF
                }

                // coarse mesh trial function loop
                for(Index j(0); j < coarse_num_loc_dofs; ++j)
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

    /// \cond internal
    namespace Intern
    {
      template<typename Matrix_>
      class RowScaler;

      template<typename DataType_>
      class RowScaler< LAFEM::SparseMatrixCSR<Mem::Main, DataType_> >
      {
      public:
        typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

        static void apply(MatrixType& matrix, const DataType_ x[])
        {
          Index* row_ptr(matrix.row_ptr());
          Index* row_end(matrix.row_ptr_end());
          DataType_* data(matrix.val());

          for(Index i(0); i < matrix.rows(); ++i)
          {
            for(Index j(row_ptr[i]); j < row_end[i]; ++j)
              data[j] *= x[i];
          }
        }
      };
    }
    /// \endcond
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_GRID_TRANSFER_HPP
