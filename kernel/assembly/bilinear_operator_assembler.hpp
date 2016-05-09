#pragma once
#ifndef KERNEL_ASSEMBLY_BILINEAR_OPERATOR_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_BILINEAR_OPERATOR_ASSEMBLER_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Bilinear Operator Assembler class template
     *
     * This class template implements the assembly of a standard bilinear operator into a matrix.
     *
     * \author Peter Zajac
     */
    class BilinearOperatorAssembler
    {
    public:
      /**
       * \brief Assembles a bilinear operator into a matrix.
       *
       * This function is the version for different test- and trial-spaces.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] operat
       * A reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] test_space
       * A reference to the finite-element test-space to be used.
       *
       * \param[in] trial_space
       * A reference to the finite-element trial-space to be used.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      template<
        typename Matrix_,
        typename Operator_,
        typename TestSpace_,
        typename TrialSpace_,
        typename CubatureFactory_>
      static void assemble_matrix2(
        Matrix_& matrix,
        Operator_& operat,
        const TestSpace_& test_space,
        const TrialSpace_& trial_space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        // matrix type
        typedef Matrix_ MatrixType;
        // functor type
        typedef Operator_ OperatorType;
        // test-space type
        typedef TestSpace_ TestSpaceType;
        // trial-space type
        typedef TrialSpace_ TrialSpaceType;

        // assembly traits
        typedef AsmTraits2<
          typename MatrixType::DataType,
          TestSpaceType,
          TrialSpaceType,
          OperatorType::trafo_config,
          OperatorType::test_config,
          OperatorType::trial_config> AsmTraits;

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

        // create a functor evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // create local matrix data
        typename AsmTraits::LocalMatrixType lmd;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // prepare functor evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_test_dofs = test_eval.get_num_local_dofs();
          int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // format local matrix
          lmd.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            test_eval(test_data, trafo_data);
            trial_eval(trial_data, trafo_data);

            // prepare bilinear operator
            oper_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_test_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_trial_dofs; ++j)
              {
                // evaluate functor and integrate
                lmd(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                  oper_eval(trial_data.phi[j], test_data.phi[i]);
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish functor evaluator
          oper_eval.finish();

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();

          // initialise dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lmd, test_dof_mapping, trial_dof_mapping, alpha);

          // finish dof mapping
          trial_dof_mapping.finish();
          test_dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }

      /**
       * \brief Assembles a bilinear operator into a matrix.
       *
       * This function is the version for identical test- and trial-spaces.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] operat
       * A reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      template<
        typename Matrix_,
        typename Operator_,
        typename Space_,
        typename CubatureFactory_>
      static void assemble_matrix1(
        Matrix_& matrix,
        Operator_& operat,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        // matrix type
        typedef Matrix_ MatrixType;
        // functor type
        typedef Operator_ OperatorType;
        // space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename MatrixType::DataType,
          SpaceType,
          OperatorType::trafo_config,
          OperatorType::test_config | OperatorType::trial_config> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functor evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local matrix data
        typename AsmTraits::LocalMatrixType lmd;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // format local matrix
          lmd.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // prepare bilinear operator
            oper_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // evaluate functor and integrate
                lmd(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                  oper_eval(space_data.phi[j], space_data.phi[i]);
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish functor evaluator
          oper_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lmd, dof_mapping, dof_mapping, alpha);

          // finish dof mapping
          dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }

      /**
       * \brief Assembles a vector valued bilinear operator into a block matrix.
       *
       * This function is the version for identical test- and trial-spaces.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled. This has to be compatible with the operator's block structure, i.e. to be
       * of some *Matrix*Blocked type with appropriate block size.
       *
       * \param[in] operat
       * A reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       *
       * \author Jordi Paul
       *
       */
      template
      <
        typename Matrix_,
        typename Operator_,
        typename Space_,
        typename CubatureFactory_
      >
      static void assemble_block_matrix1(
        Matrix_& matrix,
        Operator_& operat,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        // Make sure matrix and operator match
        static_assert(Operator_::BlockHeight == Matrix_::BlockHeight, "Operator/Matrix BlockHeight mismatch.");
        static_assert(Operator_::BlockWidth == Matrix_::BlockWidth, "Operator/Matrix BlockHeight mismatch.");
        // matrix type
        typedef Matrix_ MatrixType;
        // Block height
        static constexpr int BlockHeight = Operator_::BlockHeight;
        // Block width
        static constexpr int BlockWidth = Operator_::BlockWidth;
        // Type returned by the operator
        typedef Tiny::Matrix<typename MatrixType::DataType, BlockHeight, BlockWidth> OperatorValueType;

        // functor type
        typedef Operator_ OperatorType;
        // space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1Blocked
        <
          typename MatrixType::DataType,
          OperatorValueType,
          SpaceType,
          OperatorType::trafo_config,
          OperatorType::test_config | OperatorType::trial_config
        > AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functor evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local matrix data
        typename AsmTraits::LocalMatrixType lmd;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // format local matrix
          lmd.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // prepare bilinear operator
            oper_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // evaluate functor and integrate
                lmd(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                  oper_eval(space_data.phi[j], space_data.phi[i]);
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish functor evaluator
          oper_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lmd, dof_mapping, dof_mapping, alpha);

          // finish dof mapping
          dof_mapping.finish();

          // continue with next cell
        }
        // okay, that's it
      }

      /**
       * \brief Applies the bilinear operator to an FE function
       *
       * This function is the version for identical test- and trial-spaces and does the same as assembling the
       * operator into a matrix and multiplying it with a vector representing the FE function. Useful if this is
       * needed only a few times and one does not want to save the matrix.
       *
       * \param[out] ret
       * The return vector.
       *
       * \param[in] coeff_vector
       * The coefficient vector for the FE function.
       *
       * \param[in] operat
       * A reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       *
       * \author Jordi Paul
       *
       */
      template
      <
        typename Vector_,
        typename Operator_,
        typename Space_,
        typename CubatureFactory_
      >
      static void apply1(
        Vector_& ret,
        const Vector_& coeff_vector,
        const Operator_& operat,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Vector_::DataType alpha = typename Vector_::DataType(1))
        {
          ASSERT(ret.size()==space.get_num_dofs(), "Error: Return vector size does not match FE space dimension: " + stringify(ret.size()) + " != " + stringify(space.get_num_dofs() )) ;
          ASSERT(coeff_vector.size()==space.get_num_dofs(), "Error: Coefficient vector size does not match FE space dimension: " + stringify(coeff_vector.size()) + " != " + stringify(space.get_num_dofs() ));
          ret.format();

          // Coefficient vector type
          typedef Vector_ VectorType;
          // functor type
          typedef Operator_ OperatorType;
          // space type
          typedef Space_ SpaceType;

          // assembly traits
          typedef AsmTraits1
          <
            typename VectorType::DataType,
            SpaceType,
            OperatorType::trafo_config,
            OperatorType::test_config | OperatorType::trial_config
          > AsmTraits;

          // fetch the trafo
          const typename AsmTraits::TrafoType& trafo = space.get_trafo();

          // create a trafo evaluator
          typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

          // create a space evaluator and evaluation data
          typename AsmTraits::SpaceEvaluator space_eval(space);

          // create a dof-mapping
          typename AsmTraits::DofMapping dof_mapping(space);

          // create a functor evaluator
          typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

          // create trafo evaluation data
          typename AsmTraits::TrafoEvalData trafo_data;

          // create space evaluation data
          typename AsmTraits::SpaceEvalData space_data;

          // create local vector data
          typename AsmTraits::LocalVectorType coeff_loc;
          typename AsmTraits::LocalVectorType ret_loc;

          // create cubature rule
          typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

          // create vector scatter axpy for adding to the global return vector
          typename VectorType::ScatterAxpy scatter_axpy(ret);
          // create vector gather axpy for picking the local values from the global vector
          typename VectorType::GatherAxpy gather_axpy(coeff_vector);

          // loop over all cells of the mesh
          for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
          {
            // format local vectors
            coeff_loc.format();
            ret_loc.format();

            // initialise dof-mapping
            dof_mapping.prepare(cell);

            // prepare trafo evaluator
            trafo_eval.prepare(cell);

            // incorporate local vector
            gather_axpy(coeff_loc, dof_mapping);

            // prepare space evaluator
            space_eval.prepare(trafo_eval);

            // prepare functor evaluator
            oper_eval.prepare(trafo_eval);

            // fetch number of local dofs
            int num_loc_dofs = space_eval.get_num_local_dofs();

            // loop over all quadrature points and integrate
            for(int k(0); k < cubature_rule.get_num_points(); ++k)
            {
              // compute trafo data
              trafo_eval(trafo_data, cubature_rule.get_point(k));

              // compute basis function data
              space_eval(space_data, trafo_data);

              // prepare bilinear operator
              oper_eval.set_point(trafo_data);

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // Integrate the FE function
                  ret_loc(i) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                    oper_eval(space_data.phi[j], space_data.phi[i]) * coeff_loc(j);
                  // continue with next test function
                }
                // continue with next trial function
              }
              // continue with next cubature point
            }

            // finish functor evaluator
            oper_eval.finish();

            // finish evaluators
            space_eval.finish();
            trafo_eval.finish();

            // incorporate local vector
            scatter_axpy(ret_loc, dof_mapping, alpha);

            // finish dof mapping
            dof_mapping.finish();

            // continue with next cell
          }

          // okay, that's it
        }
      /**
       * \brief Applies the bilinear operator to an FE function
       *
       * This function is the version for different test- and trial-spaces and does the same as assembling the
       * operator into a matrix and multiplying it with a vector representing the FE function. Useful if this is
       * needed only a few times and one does not want to save the matrix.
       *
       * \param[out] ret
       * The return vector.
       *
       * \param[in] coeff_vector
       * The coefficient vector for the FE function.
       *
       * \param[in] operat
       * A reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       *
       * \author Jordi Paul
       *
       */
      template
      <
        typename Vector_,
        typename Operator_,
        typename TestSpace_,
        typename TrialSpace_,
        typename CubatureFactory_
      >
      static void apply2(
        Vector_& ret,
        const Vector_& coeff_vector,
        const Operator_& operat,
        const TestSpace_& test_space,
        const TrialSpace_& trial_space,
        const CubatureFactory_& cubature_factory,
        typename Vector_::DataType alpha = typename Vector_::DataType(1))
        {
          ASSERT(ret.size()==test_space.get_num_dofs(), "Error: Return vector size does not match FE space dimension: " + stringify(ret.size()) + " != " + stringify(test_space.get_num_dofs() ));
          ASSERT(coeff_vector.size()==trial_space.get_num_dofs(), "Error: Coefficient vector size does not match FE space dimension: " + stringify(coeff_vector.size()) + " != " + stringify(trial_space.get_num_dofs() ));
          ret.format();

          // Coefficient vector type
          typedef Vector_ VectorType;
          // functor type
          typedef Operator_ OperatorType;
          // test-space type
          typedef TestSpace_ TestSpaceType;
          // trial-space type
          typedef TrialSpace_ TrialSpaceType;

          // assembly traits
          typedef AsmTraits2
          <
            typename VectorType::DataType,
            TestSpaceType,
            TrialSpaceType,
            OperatorType::trafo_config,
            OperatorType::test_config,
            OperatorType::trial_config
          > AsmTraits;

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

          // create a functor evaluator
          typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

          // create trafo evaluation data
          typename AsmTraits::TrafoEvalData trafo_data;

          // create space evaluation data
          typename AsmTraits::TestEvalData test_data;
          typename AsmTraits::TrialEvalData trial_data;

          // create local vector data
          typename AsmTraits::LocalTrialVectorType coeff_loc;
          typename AsmTraits::LocalTestVectorType ret_loc;

          // create cubature rule
          typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

          // create vector scatter axpy for adding to the global return vector
          typename VectorType::ScatterAxpy scatter_axpy(ret);
          // create vector gather axpy for picking the local values from the global vector
          typename VectorType::GatherAxpy gather_axpy(coeff_vector);

          // loop over all cells of the mesh
          for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
          {
            // format local vectors
            coeff_loc.format();
            ret_loc.format();

            // initialise dof-mappings
            test_dof_mapping.prepare(cell);
            trial_dof_mapping.prepare(cell);

            // prepare trafo evaluator
            trafo_eval.prepare(cell);
            // prepare space evaluator
            test_eval.prepare(trafo_eval);
            trial_eval.prepare(trafo_eval);

            // incorporate local vector
            gather_axpy(coeff_loc, trial_dof_mapping);

            // prepare functor evaluator
            oper_eval.prepare(trafo_eval);

            // fetch number of local dofs
            int num_loc_test_dofs = test_eval.get_num_local_dofs();
            int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

            // loop over all quadrature points and integrate
            for(int k(0); k < cubature_rule.get_num_points(); ++k)
            {
              // compute trafo data
              trafo_eval(trafo_data, cubature_rule.get_point(k));

              // compute basis function data
              test_eval(test_data, trafo_data);
              trial_eval(trial_data, trafo_data);

              // prepare bilinear operator
              oper_eval.set_point(trafo_data);

              // test function loop
              for(int i(0); i < num_loc_test_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_trial_dofs; ++j)
                {
                  // Integrate the FE function
                  ret_loc(i) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                    oper_eval(trial_data.phi[j], test_data.phi[i]) * coeff_loc(j);
                  // continue with next trial function
                }
                // continue with next test function
              }
              // continue with next cubature point
            }

            // finish functor evaluator
            oper_eval.finish();

            // finish evaluators
            test_eval.finish();
            trial_eval.finish();
            trafo_eval.finish();

            // incorporate local vector
            scatter_axpy(ret_loc, test_dof_mapping, alpha);

            // finish dof mappings
            test_dof_mapping.finish();
            trial_dof_mapping.finish();

            // continue with next cell
          }
          // okay, that's it
        }
    }; // class BilinearOperatorAssembler<...>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_BILINEAR_OPERATOR_ASSEMBLER_HPP
