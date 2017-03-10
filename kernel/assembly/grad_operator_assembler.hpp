
#pragma once
#ifndef KERNEL_ASSEMBLY_GRAD_OPERATOR_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_GRAD_OPERATOR_ASSEMBLER_HPP 1

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Assembles the gradient operator
     *
     */
    class GradOperatorAssembler
    {
    public:

      /**
       */
      template<typename DT_, typename IT_, int dim_,
        typename SpaceVelo_, typename SpacePres_,
        typename CubatureFactory_>
      static void assemble(
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1>& matrix_g,
        const SpaceVelo_& space_velo,
        const SpacePres_& space_pres,
        const CubatureFactory_& cubature_factory,
        const DT_ scale = DT_(1)
        )
      {
        typedef DT_ DataType;
        typedef IT_ IndexType;
        static constexpr int dim = dim_;

        // ensure that the matrices have the correct dimensions
        static_assert(SpaceVelo_::shape_dim == dim_, "invalid matrix block sizes");

        // make sure that velocity and pressure have the same trafo
        if((&space_velo.get_trafo()) != (&space_pres.get_trafo()))
          throw InternalError("Velocity and Pressure spaces must be defined on the same trafo!");

        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, dim, 1> MatrixG;

        // assembly traits
        typedef AsmTraits2<
          DataType,
          SpaceVelo_,
          SpacePres_,
          TrafoTags::jac_det,
          SpaceTags::value,
          SpaceTags::grad> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space_velo.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestEvaluator velo_eval(space_velo);
        typename AsmTraits::TrialEvaluator pres_eval(space_pres);

        // create dof-mappings
        typename AsmTraits::TestDofMapping velo_dof_mapping(space_velo);
        typename AsmTraits::TrialDofMapping pres_dof_mapping(space_pres);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData velo_data;
        typename AsmTraits::TrialEvalData pres_data;

        // maximum number of dofs
        static constexpr int max_velo_dofs = AsmTraits::max_local_test_dofs;
        static constexpr int max_pres_dofs = AsmTraits::max_local_trial_dofs;

        // create local matrix data
        //typename AsmTraits::LocalMatrixType lmd;
        typedef Tiny::Matrix<DataType, dim, 1> GEntryType;
        Tiny::Matrix<GEntryType, max_velo_dofs, max_pres_dofs> local_matrix;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // first of all, check whether we need to assemble the matrix structures
        if(matrix_g.empty())
        {
          Adjacency::Graph graph_g(SymbolicAssembler::assemble_graph_std2(space_velo, space_pres));
          matrix_g = MatrixG(graph_g);
        }
        matrix_g.format();

        // create matrix scatter-axpy (if needed)
        typename MatrixG::ScatterAxpy scatter_matrix(matrix_g);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          velo_eval.prepare(trafo_eval);
          pres_eval.prepare(trafo_eval);

          // fetch number of local dofs
          const int num_loc_velo_dofs = velo_eval.get_num_local_dofs();
          const int num_loc_pres_dofs = pres_eval.get_num_local_dofs();

          // format local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int pt(0); pt < cubature_rule.get_num_points(); ++pt)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(pt));

            // compute basis function data
            velo_eval(velo_data, trafo_data);
            pres_eval(pres_data, trafo_data);

            // trial function loop
            for(int i(0); i < num_loc_velo_dofs; ++i)
            {
              // test function loop
              for(int j(0); j < num_loc_pres_dofs; ++j)
              {
                const DataType pv = trafo_data.jac_det * cubature_rule.get_weight(pt) * velo_data.phi[i].value;

                // dimension loop
                for(int k(0); k < dim; ++k)
                {
                  local_matrix[i][j][k][0] += pres_data.phi[j].grad[k] * pv;
                }

                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish evaluators
          pres_eval.finish();
          velo_eval.finish();
          trafo_eval.finish();

          // initialise dof-mappings
          velo_dof_mapping.prepare(cell);
          pres_dof_mapping.prepare(cell);

          // scatter into matrix
          scatter_g(local_matrix, velo_dof_mapping, pres_dof_mapping, scale);

          // finish dof mapping
          pres_dof_mapping.finish();
          velo_dof_mapping.finish();

          // continue with next cell
        }
      }

      template<typename DT_, typename IT_, int dim_,
        typename SpaceVelo_, typename SpacePres_,
        typename CubatureFactory_>
      static void assemble(
        LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>& vec_asm,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& vec_in,
        const SpaceVelo_& space_velo,
        const SpacePres_& space_pres,
        const CubatureFactory_& cubature_factory,
        const DT_ scale = DT_(1)
        )
      {
        typedef DT_ DataType;
        typedef IT_ IndexType;
        static constexpr int dim = dim_;

        // ensure that the matrices have the correct dimensions
        static_assert(SpaceVelo_::shape_dim == dim, "invalid matrix block sizes");

        // make sure that velocity and pressure have the same trafo
        if((&space_velo.get_trafo()) != (&space_pres.get_trafo()))
          throw InternalError("Velocity and Pressure spaces must be defined on the same trafo!");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, IndexType, dim> VectorAsm;
        typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> VectorIn;

        // assembly traits
        typedef AsmTraits2<
          DataType,
          SpaceVelo_,
          SpacePres_,
          TrafoTags::jac_det,
          SpaceTags::grad,
          SpaceTags::value> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space_velo.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestEvaluator velo_eval(space_velo);
        typename AsmTraits::TrialEvaluator pres_eval(space_pres);

        // create dof-mappings
        typename AsmTraits::TestDofMapping velo_dof_mapping(space_velo);
        typename AsmTraits::TrialDofMapping pres_dof_mapping(space_pres);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData velo_data;
        typename AsmTraits::TrialEvalData pres_data;

        // maximum number of dofs
        static constexpr int max_velo_dofs = AsmTraits::max_local_test_dofs;
        static constexpr int max_pres_dofs = AsmTraits::max_local_trial_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim> VectorValue;
        typedef Tiny::Vector<VectorValue, max_velo_dofs> LocalVectorType;
        LocalVectorType local_vector;

        typedef Tiny::Vector<DataType, max_pres_dofs> LocalVecInType;
        LocalVecInType local_vec_in_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim> loc_grad_p;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // Create vector-scatter-axpy for the vector to be assembled
        typename VectorAsm::ScatterAxpy scatter_vector(vec_asm);
        // Create vector gather axpy for the input vector
        typename VectorIn::GatherAxpy gather_vec_in(vec_in);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          velo_eval.prepare(trafo_eval);
          pres_eval.prepare(trafo_eval);

          // Initialise pressure DoF mapping
          pres_dof_mapping.prepare(cell);

          // Gather our local input DoF
          local_vec_in_dofs.format();
          gather_vec_in(local_vec_in_dofs, pres_dof_mapping);

          // fetch number of local dofs
          const int num_loc_velo_dofs = velo_eval.get_num_local_dofs();
          const int num_loc_pres_dofs = pres_eval.get_num_local_dofs();

          // format local matrix and vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int pt(0); pt < cubature_rule.get_num_points(); ++pt)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(pt));

            // compute basis function data
            velo_eval(velo_data, trafo_data);
            pres_eval(pres_data, trafo_data);

            // evaluate pressure gradient
            loc_grad_p.format();
            for(int i(0); i < num_loc_pres_dofs; ++i)
            {
              for(int k(0); k < dim; ++k)
              {
                loc_grad_p[k] += local_vec_in_dofs[i] * pres_data.phi[i].grad[k];
              }
            }

            // test function loop
            for(int i(0); i < num_loc_velo_dofs; ++i)
            {
              const DataType v(trafo_data.jac_det * cubature_rule.get_weight(pt) * velo_data.phi[i].value);

              // dimension loop
              for(int k(0); k < dim; ++k)
              {
                //local_matrix[i][j][k][0] += velo_data.phi[i].grad[k] * pv;
                // update vector
                local_vector[i][k] += v * loc_grad_p[k];
              }
            } // continue with next test function

          }// continue with next cubature point

          // finish evaluators
          pres_eval.finish();
          velo_eval.finish();
          trafo_eval.finish();

          // Initialise velocity DoF mapping
          velo_dof_mapping.prepare(cell);

          // scatter into vector
          scatter_vector(local_vector, velo_dof_mapping, scale);

          // Finish all DoF mappings
          pres_dof_mapping.finish();
          velo_dof_mapping.finish();

        } // continue with next cell
      }

    };
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_GRAD_OPERATOR_ASSEMBLER_HPP
