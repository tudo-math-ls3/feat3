// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.


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
     * This is the gradient bilinear form using blocked types without partial integration:
     * \f$ g(\psi, \phi)_m = \alpha \int_\Omega(\nabla \psi, \phi e_m) dx, m=1,\dots,d \f$
     *
     * \author Jordi Paul
     *
     */
    class GradOperatorAssembler
    {
    public:

      /**
       * \brief Assembles the bilinear form into a matrix
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam IT_
       * The index type for our containers
       *
       * \tparam dim_
       * The world dimension, which is also the block size
       *
       * \tparam SpaceTest_
       * The space for \f$ \phi \f$
       *
       * \tparam SpaceTrial_
       * The space for \f$ \psi \f$
       *
       * \tparam CubatureFactory_
       * For creating the cubature rule
       *
       * \param[out] matrix_g
       * The \transient matrix we assemble into, gets overwritten
       *
       * \param[in] test_space
       * The \transient space for \f$ \phi \f$
       *
       * \param[in] trial_space
       * The \transient space for \f$ \psi \f$
       *
       * \param[in] cubature_factory
       * Creates the cubature rule we want to use
       *
       * \param[in] scale
       * Scaling parameter \f$ \alpha \f$
       */
      template<typename DT_, typename IT_, int dim_,
        typename SpaceTest_, typename SpaceTrial_,
        typename CubatureFactory_>
      static void assemble(
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1>& matrix_g,
        const SpaceTest_& test_space,
        const SpaceTrial_& trial_space,
        const CubatureFactory_& cubature_factory,
        const DT_ scale = DT_(1)
        )
      {
        typedef DT_ DataType;
        typedef IT_ IndexType;
        static constexpr int dim = dim_;

        // ensure that the matrices have the correct dimensions
        static_assert(SpaceTest_::shape_dim == dim_, "invalid matrix block sizes");

        // make sure that trial and test spaces have the same trafo
        XASSERTM((&test_space.get_trafo()) == (&trial_space.get_trafo()),
          "Trial and test spaces must be defined on the same trafo!");

        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, dim, 1> MatrixG;

        // assembly traits
        typedef AsmTraits2<
          DataType,
          SpaceTest_,
          SpaceTrial_,
          TrafoTags::jac_det,
          SpaceTags::value,
          SpaceTags::grad> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = test_space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestEvaluator test_eval(test_space);
        typename AsmTraits::TrialEvaluator trial_eval(trial_space);

        // create dof-mappings
        typename AsmTraits::TestDofMapping test_dof_mapping(test_space);
        typename AsmTraits::TrialDofMapping trial_dof_mapping(trial_space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // maximum number of dofs
        static constexpr int max_test_dofs = AsmTraits::max_local_test_dofs;
        static constexpr int max_trial_dofs = AsmTraits::max_local_trial_dofs;

        // create local matrix data
        //typename AsmTraits::LocalMatrixType lmd;
        typedef Tiny::Matrix<DataType, dim, 1> GEntryType;
        Tiny::Matrix<GEntryType, max_test_dofs, max_trial_dofs> local_matrix;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // first of all, check whether we need to assemble the matrix structures
        if(matrix_g.empty())
        {
          Adjacency::Graph graph_g(SymbolicAssembler::assemble_graph_std2(test_space, trial_space));
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
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // fetch number of local dofs
          const int num_loc_test_dofs = test_eval.get_num_local_dofs();
          const int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // format local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int pt(0); pt < cubature_rule.get_num_points(); ++pt)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(pt));

            // compute basis function data
            test_eval(test_data, trafo_data);
            trial_eval(trial_data, trafo_data);

            // trial function loop
            for(int i(0); i < num_loc_test_dofs; ++i)
            {
              // test function loop
              for(int j(0); j < num_loc_trial_dofs; ++j)
              {
                const DataType pv = trafo_data.jac_det * cubature_rule.get_weight(pt) * test_data.phi[i].value;

                // dimension loop
                for(int k(0); k < dim; ++k)
                {
                  local_matrix[i][j][k][0] += trial_data.phi[j].grad[k] * pv;
                }

                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();

          // initialize dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // scatter into matrix
          scatter_matrix(local_matrix, test_dof_mapping, trial_dof_mapping, scale);

          // finish dof mapping
          trial_dof_mapping.finish();
          test_dof_mapping.finish();

          // continue with next cell
        }
      }

      /**
       * \brief Assembles the application of the bilinear form to a vector into a vector
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam IT_
       * The index type for our containers
       *
       * \tparam dim_
       * The world dimension, which is also the block size
       *
       * \tparam SpaceTest_
       * The space for \f$ \phi \f$
       *
       * \tparam SpaceTrial_
       * The space for \f$ \psi \f$
       *
       * \tparam CubatureFactory_
       * For creating the cubature rule
       *
       * \param[in, out] vec_asm
       * The \transient vector we assemble onto
       *
       * \param[in] vec_in
       * The \transient vector we apply the bilinear form to
       *
       * \param[in] test_space
       * The \transient space for \f$ \phi \f$
       *
       * \param[in] trial_space
       * The \transient space for \f$ \psi \f$
       *
       * \param[in] cubature_factory
       * Creates the cubature rule we want to use
       *
       * \param[in] scale
       * Scaling parameter \f$ \alpha \f$
       */
      template<typename DT_, typename IT_, int dim_,
        typename SpaceTest_, typename SpaceTrial_,
        typename CubatureFactory_>
      static void assemble(
        LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>& vec_asm,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& vec_in,
        const SpaceTest_& test_space,
        const SpaceTrial_& trial_space,
        const CubatureFactory_& cubature_factory,
        const DT_ scale = DT_(1)
        )
      {
        typedef DT_ DataType;
        typedef IT_ IndexType;
        static constexpr int dim = dim_;

        // ensure that the matrices have the correct dimensions
        static_assert(SpaceTest_::shape_dim == dim, "invalid matrix block sizes");

        // make sure that trial and test spaces have the same trafo
        XASSERTM((&test_space.get_trafo()) == (&trial_space.get_trafo()),
          "Trial and test spaces must be defined on the same trafo!");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType, IndexType, dim> VectorAsm;
        typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> VectorIn;

        // assembly traits
        typedef AsmTraits2<
          DataType,
          SpaceTest_,
          SpaceTrial_,
          TrafoTags::jac_det,
          SpaceTags::value,
          SpaceTags::grad> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = test_space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestEvaluator test_eval(test_space);
        typename AsmTraits::TrialEvaluator trial_eval(trial_space);

        // create dof-mappings
        typename AsmTraits::TestDofMapping test_dof_mapping(test_space);
        typename AsmTraits::TrialDofMapping trial_dof_mapping(trial_space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // maximum number of dofs
        static constexpr int max_test_dofs = AsmTraits::max_local_test_dofs;
        static constexpr int max_trial_dofs = AsmTraits::max_local_trial_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim> VectorValue;
        typedef Tiny::Vector<VectorValue, max_test_dofs> LocalVectorType;
        LocalVectorType local_vector;

        typedef Tiny::Vector<DataType, max_trial_dofs> LocalVecInType;
        LocalVecInType local_vec_in_dofs;

        // our local trial function gradient
        Tiny::Vector<DataType, dim> loc_grad_trial;

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
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // Initialize trial DoF mapping
          trial_dof_mapping.prepare(cell);

          // Gather our local input DoF
          local_vec_in_dofs.format();
          gather_vec_in(local_vec_in_dofs, trial_dof_mapping);

          // fetch number of local dofs
          const int num_loc_test_dofs = test_eval.get_num_local_dofs();
          const int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // format local matrix and vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int pt(0); pt < cubature_rule.get_num_points(); ++pt)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(pt));

            // compute basis function data
            test_eval(test_data, trafo_data);
            trial_eval(trial_data, trafo_data);

            // evaluate trial function gradient
            loc_grad_trial.format();
            for(int i(0); i < num_loc_trial_dofs; ++i)
            {
              for(int k(0); k < dim; ++k)
              {
                loc_grad_trial[k] += local_vec_in_dofs[i] * trial_data.phi[i].grad[k];
              }
            }

            // test function loop
            for(int i(0); i < num_loc_test_dofs; ++i)
            {
              const DataType v(trafo_data.jac_det * cubature_rule.get_weight(pt) * test_data.phi[i].value);

              // dimension loop
              for(int k(0); k < dim; ++k)
              {
                //local_matrix[i][j][k][0] += test_data.phi[i].grad[k] * pv;
                // update vector
                local_vector[i][k] += v * loc_grad_trial[k];
              }
            } // continue with next test function

          }// continue with next cubature point

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();

          // Initialize test DoF mapping
          test_dof_mapping.prepare(cell);

          // scatter into vector
          scatter_vector(local_vector, test_dof_mapping, scale);

          // Finish all DoF mappings
          trial_dof_mapping.finish();
          test_dof_mapping.finish();

        } // continue with next cell
      }

    };
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_GRAD_OPERATOR_ASSEMBLER_HPP
