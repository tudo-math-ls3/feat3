// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/likwid_marker.hpp>

namespace FEAT
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
       * The \transient matrix that is to be assembled.
       *
       * \param[in] operat
       * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] test_space
       * A \transient reference to the finite-element test-space to be used.
       *
       * \param[in] trial_space
       * A \transient reference to the finite-element trial-space to be used.
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory to be used for integration.
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
        // validate matrix dimensions
        XASSERTM(matrix.rows() == test_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == trial_space.get_num_dofs(), "invalid matrix dimensions");

        // matrix type
        typedef Matrix_ MatrixType;
        // operator type
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

        // create a operator evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // the value type of the operator
        typedef typename OperatorType::template Evaluator<AsmTraits>::ValueType ValueType;

        // ensure that the operator and matrix value types are compatible
        static_assert(std::is_same<ValueType, typename MatrixType::ValueType>::value,
          "matrix and bilinear operator have different value types!");

        // create local matrix data
        typename AsmTraits::template TLocalMatrix<ValueType> loc_mat;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);
        FEAT_APPLICATION_MARKER_START("bilin_op_asm_mat2");
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // prepare operator evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_test_dofs = test_eval.get_num_local_dofs();
          int num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // format local matrix
          loc_mat.format();

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
                // evaluate operator and integrate
                Tiny::axpy(loc_mat(i,j), oper_eval.eval(trial_data.phi[j], test_data.phi[i]),
                  trafo_data.jac_det * cubature_rule.get_weight(k));
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }
          FEAT_APPLICATION_MARKER_STOP("bilin_op_asm_mat2");

          // finish operator evaluator
          oper_eval.finish();

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();

          // initialize dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(loc_mat, test_dof_mapping, trial_dof_mapping, alpha);

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
       * The \transient matrix that is to be assembled.
       *
       * \param[in] operat
       * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A \transient reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory to be used for integration.
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
        // validate matrix dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");

        // matrix type
        typedef Matrix_ MatrixType;
        // operator type
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

        // create a operator evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // the value type of the operator
        typedef typename OperatorType::template Evaluator<AsmTraits>::ValueType ValueType;

        // ensure that the operator and matrix value types are compatible
        static_assert(std::is_same<ValueType, typename MatrixType::ValueType>::value,
          "matrix and bilinear operator have different value types!");

        // create local matrix data
        typename AsmTraits::template TLocalMatrix<ValueType> loc_mat;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_axpy(matrix);
        FEAT_APPLICATION_MARKER_START("bilin_op_asm:matrix1");

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare operator evaluator
          oper_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // format local matrix
          loc_mat.format();

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
                // evaluate operator and integrate
                Tiny::axpy(loc_mat(i,j), oper_eval.eval(space_data.phi[j], space_data.phi[i]),
                  trafo_data.jac_det * cubature_rule.get_weight(k));
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }
        FEAT_APPLICATION_MARKER_STOP("bilin_op_asm:matrix1");
          // finish operator evaluator
          oper_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(loc_mat, dof_mapping, dof_mapping, alpha);

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
       * The \transient return vector.
       *
       * \param[in] coeff_vector
       * The \transient coefficient vector for the FE function.
       *
       * \param[in] operat
       * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] space
       * A \transient reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       *
       * \author Jordi Paul
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
        XASSERTM(ret.size()==space.get_num_dofs(), "Return vector size does not match FE space dimension");
        XASSERTM(coeff_vector.size()==space.get_num_dofs(), "Coefficient vector size does not match FE space dimension");
        ret.format();

        // Coefficient vector type
        typedef Vector_ VectorType;
        // operator type
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

        // create a operator evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // get the value type from the vector
        typedef typename VectorType::ValueType ValueType;

        // create local vector data
        typename AsmTraits::template TLocalVector<ValueType> loc_coeff, loc_ret;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector scatter axpy for adding to the global return vector
        typename VectorType::ScatterAxpy scatter_axpy(ret);
        // create vector gather axpy for picking the local values from the global vector
        typename VectorType::GatherAxpy gather_axpy(coeff_vector);
        FEAT_APPLICATION_MARKER_START("bilin_op_asm:apply1");
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // format local vectors
          loc_coeff.format();
          loc_ret.format();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // incorporate local vector
          gather_axpy(loc_coeff, dof_mapping);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare operator evaluator
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
                Tiny::axpy(loc_ret(i),
                  oper_eval.eval(space_data.phi[j], space_data.phi[i]) * loc_coeff(j),
                  trafo_data.jac_det * cubature_rule.get_weight(k));
                // continue with next test function
              }
              // continue with next trial function
            }
            // continue with next cubature point
          }
        FEAT_APPLICATION_MARKER_STOP("bilin_op_asm:apply1");
          // finish operator evaluator
          oper_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // incorporate local vector
          scatter_axpy(loc_ret, dof_mapping, alpha);

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
       * The \transient return vector.
       *
       * \param[in] coeff_vector
       * The \transient coefficient vector for the FE function.
       *
       * \param[in] operat
       * A \transient reference to the operator implementing the BilinearOperator interface to be assembled.
       *
       * \param[in] test_space
       * A \transient reference to the finite-element test-space to be used.
       *
       * \param[in] trial_space
       * A \transient reference to the finite-element trial-space to be used.
       *
       * \param[in] cubature_factory
       * A \transient reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       *
       * \author Jordi Paul
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
        XASSERTM(ret.size()==test_space.get_num_dofs(), "Return vector size does not match FE space dimension");
        XASSERTM(coeff_vector.size()==trial_space.get_num_dofs(), "Coefficient vector size does not match FE space dimension");
        ret.format();

        // Coefficient vector type
        typedef Vector_ VectorType;
        // operator type
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

        // create a operator evaluator
        typename OperatorType::template Evaluator<AsmTraits> oper_eval(operat);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;
        typename AsmTraits::TrialEvalData trial_data;

        // get the value type from the vector
        typedef typename VectorType::ValueType ValueType;

        // create local vector data
        typename AsmTraits::template TLocalTrialVector<ValueType> loc_coeff;
        typename AsmTraits::template TLocalTestVector<ValueType> loc_ret;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector scatter axpy for adding to the global return vector
        typename VectorType::ScatterAxpy scatter_axpy(ret);
        // create vector gather axpy for picking the local values from the global vector
        typename VectorType::GatherAxpy gather_axpy(coeff_vector);

        FEAT_APPLICATION_MARKER_START("bilin_op_asm:apply2");
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // format local vectors
          loc_coeff.format();
          loc_ret.format();

          // initialize dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);
          // prepare space evaluator
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // incorporate local vector
          gather_axpy(loc_coeff, trial_dof_mapping);

          // prepare operator evaluator
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
                Tiny::axpy(loc_ret(i),
                  oper_eval.eval(trial_data.phi[j], test_data.phi[i]) * loc_coeff(j),
                  trafo_data.jac_det * cubature_rule.get_weight(k));
                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }
        FEAT_KERNEL_MARKER_STOP("bilin_op_asm:apply2");

          // finish operator evaluator
          oper_eval.finish();

          // finish evaluators
          test_eval.finish();
          trial_eval.finish();
          trafo_eval.finish();

          // incorporate local vector
          scatter_axpy(loc_ret, test_dof_mapping, alpha);

          // finish dof mappings
          test_dof_mapping.finish();
          trial_dof_mapping.finish();

          // continue with next cell
        }
        // okay, that's it
      }
    }; // class BilinearOperatorAssembler<...>
  } // namespace Assembly
} // namespace FEAT
